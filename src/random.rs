use std::{sync::Arc, collections::HashMap, path::PathBuf};

use rand::{thread_rng, seq::SliceRandom};
use rayon::prelude::{IntoParallelRefIterator, ParallelIterator};

use crate::{contigs::Contig, bin_classes::BinGenerator, evo_alg_utils::generate_random_permutation, bin_info_storing::{BinSet, BinInfoStorer, BinScorer}, utils::test_permutation_himem};



fn random_permutation_test(contig_sets: Vec<Vec<Arc<Contig>>>, number_of_random_permutations: usize, bin_generator: Arc<BinGenerator>, bin_scorer: Arc<BinScorer>) {

    let contig_set_id_hashmap = create_hashmap_of_contig_id_to_contig_sets(contig_sets);
    let vec_of_bin_ids: Vec<usize> = contig_set_id_hashmap.keys().cloned().collect();

    let mut all_random_permuts_to_test = Vec::new();

    for i in 0..number_of_random_permutations {
        all_random_permuts_to_test.push(generate_random_permutation(&vec_of_bin_ids, &all_random_permuts_to_test));
    }


}


fn create_hashmap_of_contig_id_to_contig_sets(contig_sets: Vec<Vec<Arc<Contig>>>) -> HashMap<usize, Vec<Arc<Contig>>> {
    let mut contig_set_id_hashmap = HashMap::new();
    for (current_id, contig_set) in contig_sets.into_iter().enumerate() {
        contig_set_id_hashmap.insert(current_id, contig_set);
    }
    contig_set_id_hashmap
}

fn construct_randomised_binset(contigs: &mut Vec<Arc<Contig>>) -> Vec<Vec<Arc<Contig>>> {
    let mut rng = thread_rng();

    contigs.shuffle(&mut rng);
    let mut randomised_bin_set = Vec::new();
    let mut current_bin = Vec::new();
    
    for contig in contigs {
        
        current_bin.push(contig.clone());
        
        if rand::random() {
            randomised_bin_set.push(current_bin);
            current_bin = Vec::new();
        }
    }
    randomised_bin_set.push(current_bin);
    randomised_bin_set

}


fn test_the_permutations(bin_generator: Arc<BinGenerator>, contig_set_id_hashmap: HashMap<usize, Vec<Arc<Contig>>>, permutations_by_id: Vec<Vec<usize>>, threads: usize, bin_scorer: Arc<BinScorer>) -> Vec<BinSet> {
    let mut permutations_by_contigs_sets = convert_permutations_by_id_to_permutations_by_contig_set(permutations_by_id, contig_set_id_hashmap).into_iter();

    let bin_info_storer = BinInfoStorer::initialise_bin_info_storer();
    let arc_bin_info_storer = Arc::new(bin_info_storer);
    let mut tested_all_permutations = false;
    let mut all_generated_bin_sets: Vec<BinSet> = Vec::new();
    while tested_all_permutations == false {

        let mut round_of_permutations = Vec::new();

        for i in 0..threads {
            match permutations_by_contigs_sets.next() {

                Some(permut_to_test) => round_of_permutations.push(permut_to_test),
                None => tested_all_permutations = true

            };
        }
        let bin_sets_generated: Vec<BinSet> = round_of_permutations.par_iter().map(|permutation| {
            let arc_bin_scorer_clone = Arc::clone(&bin_scorer);
            test_permutation_himem(permutation.clone(), Arc::clone(&bin_generator), Arc::clone(&arc_bin_info_storer), arc_bin_scorer_clone)
        }).collect();
        all_generated_bin_sets.extend(bin_sets_generated);
    };

    all_generated_bin_sets


}


fn convert_permutations_by_id_to_permutations_by_contig_set(permutations_by_id: Vec<Vec<usize>>, contig_set_id_hashmap: HashMap<usize, Vec<Arc<Contig>>>) -> Vec<Vec<Vec<Arc<Contig>>>> {
    let mut permutations_by_contig_sets = Vec::new();
    for permutation in permutations_by_id {
        let permutation_by_contig_set: Vec<Vec<Arc<Contig>>> = permutation.iter().map(|id| contig_set_id_hashmap.get(id).unwrap().clone()).collect();
        permutations_by_contig_sets.push(permutation_by_contig_set);
    }
    permutations_by_contig_sets
}

fn generate_result_information(hash_directory: &PathBuf, results_directory: &PathBuf, generated_bin_sets: Vec<BinSet>) {
    for (bin_num, binset) in generated_bin_sets.iter().enumerate() {
        let best_bin_directory = results_directory.join(&format!("{}", bin_num));
        binset.create_best_bin_dir_and_info_from_best_hashes(hash_directory, &best_bin_directory, false);
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use std::{fs, env, path::Path};
    use lazy_static::lazy_static;

    #[test]
    fn run_full_random_permutation_test() {

    }

}
