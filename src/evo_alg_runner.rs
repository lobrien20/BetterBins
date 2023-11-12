use std::{sync::Arc, collections::HashMap, path::PathBuf, fs::File, io::Write, cmp::Ordering};

use itertools::enumerate;
use log::{debug, info};
use non_dominated_sort::DominanceOrd;
use rand::{Rng, seq::SliceRandom};
use rayon::prelude::{IntoParallelRefIterator, ParallelIterator};

use crate::{bin_info_storing::{BinSet, BinScorer, BinInfoStorer}, contigs::Contig, bin_classes::BinGenerator, evo_alg_gen_tracking::EvoAlgRunner, evo_alg_utils::{mutate_permutation, get_parent_picking_approach_from_string, OffSpringCreator, ChildGenerator}, utils::test_permutation_himem};
use non_dominated_sort::non_dominated_sort;

pub fn initialise_evolutionary_alg_structs(perc_mutation: f64, initial_potential_bins: Vec<Vec<Arc<Contig>>>, perc_crossover: f64, permutations_per_generation: usize, max_generations: usize, parent_picker: String) -> EvoAlgRunner {
    info!("Using {} total bins", initial_potential_bins.len());

    let child_generator = ChildGenerator::initialise(perc_mutation, initial_potential_bins.len(), perc_crossover);
    match &child_generator.crossover {
        Some(_) => {
            if permutations_per_generation % 2 != 0 {
                panic!("Using parent crossover with odd number of permutations. Pick an even number")
            }

        },
        None => ()
    }
    let parent_picker = get_parent_picking_approach_from_string(parent_picker, permutations_per_generation);
    let evo_alg_runner = EvoAlgRunner::initialise(initial_potential_bins, parent_picker, permutations_per_generation, max_generations, child_generator);
    info!("Total generations is {}, permutations per generation is: {}", evo_alg_runner.max_generations, evo_alg_runner.permutations_per_generation);

    evo_alg_runner
}

pub fn run_evolutionary_alg(bin_generator: Arc<BinGenerator>, bin_scorer: Arc<BinScorer>, res_dir: &PathBuf, hash_directory: &PathBuf, mut evo_alg_runner: EvoAlgRunner) {
    let mut bin_info_storer = BinInfoStorer::initialise_bin_info_storer();

    let mut new_permutations_by_id = evo_alg_runner.create_initial_permutations();
    let bin_info_store_ref = Arc::new(bin_info_storer);
    let mut completed_bin_sets = test_generation_himem(new_permutations_by_id, bin_generator.clone(), bin_scorer.clone(), &mut evo_alg_runner, Arc::clone(&bin_info_store_ref));
    bin_info_storer = Arc::try_unwrap(bin_info_store_ref).unwrap();
    bin_info_storer.add_bins_from_bin_sets_to_arc_contig_hashmap(&completed_bin_sets);
    while evo_alg_runner.evo_alg_state.generation != evo_alg_runner.max_generations {

        debug!("Entering generation: {}", &evo_alg_runner.evo_alg_state.generation);
        new_permutations_by_id = evo_alg_runner.create_new_generation(completed_bin_sets);
        let bin_info_store_ref = Arc::new(bin_info_storer);
        completed_bin_sets = test_generation_himem(new_permutations_by_id, bin_generator.clone(), bin_scorer.clone(), &mut evo_alg_runner, Arc::clone(&bin_info_store_ref));
        bin_info_storer = Arc::try_unwrap(bin_info_store_ref).unwrap();
        bin_info_storer.add_bins_from_bin_sets_to_arc_contig_hashmap(&completed_bin_sets);


    }

    process_results(evo_alg_runner, bin_scorer, res_dir, hash_directory);

}






fn test_generation_himem(permutations_by_bin_set_id: Vec<Vec<usize>>, bin_generator: Arc<BinGenerator>, bin_scorer: Arc<BinScorer>, evo_alg_runner: &mut EvoAlgRunner, bin_info_storer: Arc<BinInfoStorer>) -> Vec<BinSet> {
    debug!("Testing generation! {} permutations to test.", permutations_by_bin_set_id.len());
    let permutations_to_test: Vec<(Vec<usize>, Vec<Vec<Arc<Contig>>>)> = permutations_by_bin_set_id.iter().map(|permutation| (permutation.clone(), evo_alg_runner.generate_contig_set_order_from_permutation(permutation))).collect();
    debug!("{} unique permutations to test (mapped)", permutations_to_test.len());
    let completed_bin_sets: Vec<BinSet> = permutations_to_test.par_iter().map(|(bin_order, permutation)| {
        let bin_info_store_ref = Arc::clone(&bin_info_storer);
        let bin_gen_ref = Arc::clone(&bin_generator);
        let bin_score_ref = Arc::clone(&bin_scorer);
        let mut bin_set = test_permutation_himem(permutation.clone(), bin_gen_ref, bin_info_store_ref, bin_score_ref);

        bin_set.bin_set_order = Some(bin_order.clone());
        bin_set.overall_score = Some(bin_scorer.score_bin_set(&bin_set));
        bin_set.generate_overall_comp_and_cont(&bin_scorer);
        bin_set


    }).collect();



    info!("Generation test complete. Successfully tested {} permutations.", completed_bin_sets.len());

    evo_alg_runner.evo_alg_state.next_generation(&completed_bin_sets);
    completed_bin_sets
}







fn process_results(evo_alg_runner: EvoAlgRunner, bin_scorer: Arc<BinScorer>, res_dir: &PathBuf, hash_directory: &PathBuf) {

    let mut bin_sets_string = "bin_set_number\tbin_set_score\toverall_completeness\toverall_contamination\toverall_completeness(weighted)\toverall_contamination(weighted)\n".to_string();
    let mut best_bin_sets = bin_scorer.get_best_bin_sets(&evo_alg_runner.evo_alg_state.used_permutations);
    let (least_contam, most_complete) = bin_scorer.get_least_contam_and_most_complete(&best_bin_sets);
    let mut bin_set_num_hashmap = HashMap::new();
    for (permut_num, bin_set) in enumerate(evo_alg_runner.evo_alg_state.used_permutations) {
        let mut permutation_res_dir = res_dir.join(&format!("{}_results/", permut_num));
        if best_bin_sets.contains(&bin_set) {
            bin_set_num_hashmap.insert(permut_num, bin_set.clone());
            let mut permutation_res_dir = res_dir.join(&format!("{}_results/", permut_num));
            bin_set.create_best_bin_dir_and_info_from_best_hashes(hash_directory, &permutation_res_dir, true);

        }

        bin_sets_string.push_str(&format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\n", permut_num, bin_set.overall_score.unwrap(), bin_set.best_bins.len(), bin_set.get_unweighted_completion(), bin_set.get_unweighted_contamination(), bin_set.overall_completeness, bin_set.overall_contamination).to_string());


    }


    let mut binset_info_path = res_dir.join("binset_info_path.tsv");
    let mut bin_info_file = File::create(binset_info_path).unwrap();
    bin_info_file.write_all(bin_sets_string.as_bytes());

}

















mod tests {
    use super::*;
    use std::{fs, env, path::Path};

    fn create_some_test_permutations() -> Vec<Vec<usize>> {

        let permut_1 = vec![0, 1, 2, 3, 4, 6, 5];
        let permut_2 = vec![2, 4, 3, 1, 5, 0, 6];
        let permut_3 = vec![6, 4, 1, 5, 3, 0, 2];
        let permut_4 = vec![1, 4, 2, 5, 0, 3, 6];

        vec![permut_1, permut_2, permut_3, permut_4]
    }
    fn create_some_test_best_bins() -> Vec<BinSet> {
        let test_permutations = create_some_test_permutations();
        let mut test_best_bins = Vec::new();
        let mut score_to_add = 0.0;
        for permut in test_permutations {

            let mut best_bins_t= BinSet::constructnew();

            best_bins_t.overall_score = Some(score_to_add);
            best_bins_t.bin_set_order = Some(permut);
            test_best_bins.push(best_bins_t);
            score_to_add += 1.0;

        }
        test_best_bins
    }
   /*  fn create_test_evo_alg_state() -> EvoAlgState {
        let test_best_bins = create_some_test_best_bins();

        EvoAlgState { generation: 0, used_permutations: test_best_bins.clone(), best_permutation: Some(test_best_bins[test_best_bins.len() - 1].clone()) }
    } */
  /*  fn create_test_evo_alg_runner() -> EvoAlgRunner {
        let parent_picker = get_parent_picking_approach_from_string("Elite".to_string(), 40);
        EvoAlgRunner::initialise(Vec::new(), parent_picker)


    } */
    fn create_fake_contig_sets() -> Vec<Vec<Arc<Contig>>> {

    let fake_contig_1 = Contig::new_contig(">fake_contig_1".to_string(), "ATGGCTAGCATCGATGCTAGCAGGAGCGGAGAGCTATGCATGC\n".to_string());
    let fake_contig_2 = Contig::new_contig(">fake_contig_2".to_string(), "ATGGCTAGCATCGATGCTTTTTTAGAGCTATGCATGC\n".to_string());
    let fake_contig_3 = Contig::new_contig(">fake_contig_3".to_string(), "AGTTTTTTGCTAGCATCGATGCTTTTTTGGGGGGC\n".to_string());
    let fake_contig_4 = Contig::new_contig(">fake_contig_4".to_string(), "ATGGCTAGCATCGATGCTTTTTTAGAGCTATGCATGC\n".to_string());
    let fake_contig_5 = Contig::new_contig(">fake_contig_5".to_string(), "ATGGCTAGCAAAATCGATGCTTTTTTAGAGCTATGCATGC\n".to_string());
    let fake_contig_6 = Contig::new_contig(">fake_contig_6".to_string(), "ATGGATATATATCTAGCATCGATGCTTTGTTAGAGCTACCCTGCATGC\n".to_string());
    let fake_contig_7 = Contig::new_contig(">fake_contig_6".to_string(), "ATGGCTACCCGCATCGATGCTTTGTTAGAGCTACCCTGCATGC\n".to_string());
    let fake_contig_8 = Contig::new_contig(">fake_contig_6".to_string(), "ATGGGGGGCTAGCATCGATGCTTTGTTAGAGCTACCCTGCATGC\n".to_string());

    let fake_contig_set_1 = vec![Arc::new(fake_contig_1)];
    let fake_contig_set_2 = vec![Arc::new(fake_contig_2)];
    let fake_contig_set_3 = vec![Arc::new(fake_contig_3)];
    let fake_contig_set_4 = vec![Arc::new(fake_contig_4)];
    let fake_contig_set_5 = vec![Arc::new(fake_contig_5)];
    let fake_contig_set_6 = vec![Arc::new(fake_contig_6)];
    let fake_contig_set_7 = vec![Arc::new(fake_contig_7), Arc::new(fake_contig_8)];

    vec![fake_contig_set_1, fake_contig_set_2, fake_contig_set_3, fake_contig_set_4, fake_contig_set_5, fake_contig_set_6, fake_contig_set_7]
    }

/*
    #[test]
    fn test_create_new_generation() {

        let mut evo_alg_runner = create_test_evo_alg_runner();
        let test_best_bins = create_some_test_best_bins();
        let test_generation_result = evo_alg_runner.create_new_generation(create_some_test_best_bins());
        assert_eq!(test_generation_result.len(), 2);
        let best_bin_order = test_best_bins[test_best_bins.len() - 1].bin_set_order.clone().unwrap();
        let second_best_bin_order = test_best_bins[test_best_bins.len() - 1].bin_set_order.clone().unwrap();
        let mut first_bin_closeness = 0;
        let mut second_bin_closeness = 0;
        for ((test_b, first_bin), second_bin) in test_generation_result[0].iter().zip(best_bin_order.iter()).zip(second_best_bin_order.iter()) {
            if test_b == first_bin {
                first_bin_closeness += 1;
            }
            if test_b == second_bin {
                second_bin_closeness += 1;
            }
        }
        println!("{:?}", &best_bin_order);
        println!("{:?}", &test_generation_result[0]);
        assert_eq!(first_bin_closeness, (best_bin_order.len() - 2));
        assert_eq!(second_bin_closeness, (best_bin_order.len() - 2));


    }

*/

    /*  fn test_create_bin_sets_from_original() {
        let parent_picker = get_parent_picking_approach_from_string("Elite".to_string(), 40);
        let mut evo_alg_runner = EvoAlgRunner::initialise(create_fake_contig_sets(), parent_picker);


        let mut test_original = create_some_test_best_bins()[0].clone();
        test_original.bin_set_order = Some(vec![7, 1, 2]);
    //    println!("{:?}", test_original.)
        println!("{:?}", test_original.clone().bin_set_order.unwrap());
        let new_permutations = evo_alg_runner.create_new_bin_set_permutations_from_original(2, &mut test_original);
        println!("{:?}", new_permutations);

    }
    */

}
