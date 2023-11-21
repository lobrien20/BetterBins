use std::collections::{HashMap, HashSet};
use std::path::PathBuf;
use std::sync::Arc;

use itertools::Itertools;
use log::info;
use rayon::prelude::{IntoParallelRefIterator, IntoParallelIterator};
use rayon::iter::ParallelIterator;
use crate::bin_generator::BinGenerator;
use crate::bin_info_storage::Bin;
use crate::{contigs::Contig, bin_generator::BinGen};
use crate::{bin_sets::BinSet};
use crate::bin_scoring::BinScorer;




pub fn run_classic_best_bin(mut contig_sets: Vec<Vec<Arc<Contig>>>, bin_generator: Box<dyn BinGenerator>, bin_scorer: &BinScorer) -> Vec<Bin> {
    let mut current_best_bins = Vec::new();
    if contig_sets.len() == 0 {
        panic!("No bins to run classic best bin on.");
    }
    loop {
        match generate_hypothetical_bins_from_contig_sets(contig_sets, &bin_generator) {
            
            Some(potential_bins) => { 
                
                let best_bin = identify_best_bin_from_potential_bin_choices(&potential_bins, bin_scorer);
                let best_bin_contigs = best_bin.bin_contigs.clone();
                current_best_bins.push(best_bin);
                
                match generate_new_contig_sets_with_current_best_best_bin_removed_from_current_bins(potential_bins, best_bin_contigs) {
                    
                    Some(current_cycle_contig_set) => contig_sets = current_cycle_contig_set,
                    None => break
                
                    }
                
                },
            
            None => break
            
            }

        }
    
    info!("No more remaining bins!");
    current_best_bins


    }

fn generate_hypothetical_bins_from_contig_sets(contig_sets: Vec<Vec<Arc<Contig>>>, bin_generator: &Box<dyn BinGenerator>) -> Option<Vec<Bin>> {
    let generated_bins: Vec<Bin> = contig_sets.into_par_iter().filter_map(|contig_set| {
        match bin_generator.generate_new_bin_from_contigs(contig_set.clone()) {
            Some(bin) => Some(bin),
            None => None,
        }
         }).collect();

    if generated_bins.is_empty() {
        None
    } else {
        Some(generated_bins)
    }
}
    
fn generate_new_contig_sets_with_current_best_best_bin_removed_from_current_bins(bins: Vec<Bin>, best_bin_contigs: Vec<Arc<Contig>>) -> Option<Vec<Vec<Arc<Contig>>>> {
    
    let filtered_contig_sets = bins.into_iter()
        .map(|bin| {
            bin.bin_contigs.into_iter()
            .filter(|contig| !best_bin_contigs.contains(&contig)).collect_vec()}
        ).filter(|contig_set| contig_set.len() > 0).collect_vec();
    
    info!("{} remaining potential bins", filtered_contig_sets.len());
    if filtered_contig_sets.len() > 0 {

        Some(filtered_contig_sets)
    } else {
        None
    }
}


fn identify_best_bin_from_potential_bin_choices(potential_bins: &Vec<Bin>, bin_scorer: &BinScorer) -> Bin {
    
    let mut bin_score_vec: Vec<(&Bin, f64)> = potential_bins.into_iter()
        .map(|bin| (bin, bin_scorer.score_bin(bin.completeness, bin.contamination)))
        .collect_vec();
    
    bin_score_vec.sort_by(|(_, a), (_, b)| b.partial_cmp(a).unwrap());
    


    let mut potential_best_bins = bin_score_vec.iter()
        .filter(|(_, bin_score)| bin_score == &bin_score_vec[0].1).collect_vec();

        if potential_best_bins.len() > 1 {
            potential_best_bins.sort_by(|a, b| a.0.bin_contigs.len().cmp(&b.0.bin_contigs.len()));
            info!("Multiple best bins of same score, bin picked instead by least amount of contigs!");

    }
    info!("Best bin has score of {} and {} contigs", potential_best_bins[0].1, potential_best_bins[0].0.bin_contigs.len());
    potential_best_bins[0].0.clone()


}






mod tests {
    use crate::{bin_info_storage::BinType, bin_generator};

    use super::*;
    use std::{fs, env, path::Path};
    use itertools::Itertools;
    use lazy_static::lazy_static;

    lazy_static! {
        static ref TEST_MODULE_DIR: PathBuf = PathBuf::from("tests/new_tests/classic_best_bin_test/");
        static ref CHECKM2_DB_PATH: PathBuf = PathBuf::from("tests/new_tests/databases_for_testing/uniref100.KO.1.dmnd");
        static ref COMPLEASM_DB_LIB: PathBuf = PathBuf::from("tests/new_tests/databases_for_testing/eukaryota_odb10");
    }


    fn create_fake_contigs_for_unit_test() -> Vec<Arc<Contig>> {
        let fake_contig_1 = Contig::new_contig(">fake_contig_1".to_string(), "ATGGCTAGCATCGATGCTAGCAGGAGCGGAGAGCTATGCATGC\n".to_string());
        let fake_contig_2 = Contig::new_contig(">fake_contig_2".to_string(), "ATGGCTAGCATCGATGCTTTTTTAGAGCTATGCATGC\n".to_string());
        let fake_contig_3 = Contig::new_contig(">fake_contig_3".to_string(), "AGTTTTTTGCTAGCATCGATGCTTTTTTGGGGGGC\n".to_string());
        let fake_contig_4 = Contig::new_contig(">fake_contig_4".to_string(), "ATGGCTAGCATCGATGCTTTTTTAGAGCTATGCATGC\n".to_string());
        let fake_contig_5 = Contig::new_contig(">fake_contig_5".to_string(), "ATGGCTAGCAAAATCGATGCTTTTTTAGAGCTATGCATGC\n".to_string());
        let fake_contig_6 = Contig::new_contig(">fake_contig_6".to_string(), "ATGGATATATATCTAGCATCGATGCTTTGTTAGAGCTACCCTGCATGC\n".to_string());
        let fake_contig_7 = Contig::new_contig(">fake_contig_7".to_string(), "ATGGCTACCCGCATCGATGCTTTGTTAGAGCTACCCTGCATGC\n".to_string());
        let fake_contig_8 = Contig::new_contig(">fake_contig_8".to_string(), "ATGGGGGGCTAGCATCGATGCTTTGTTAGAGCTACCCTGCATGC\n".to_string());
        let contig_vec = vec![fake_contig_1, fake_contig_2, fake_contig_3, fake_contig_4, fake_contig_5, fake_contig_6, fake_contig_7, fake_contig_8];
        let fake_contig_arc: Vec<Arc<Contig>> = contig_vec.into_iter().map(|x| Arc::new(x)).collect();
        fake_contig_arc
    }
    /* 
    fn create_fake_contig_set_bins() -> Vec<Bin> {
        let fake_contigs = create_fake_contigs_for_unit_test();
        let fake_first_best_contigs = vec![Arc::clone(&fake_contigs[0]), Arc::clone(&fake_contigs[1]), Arc::clone(&fake_contigs[2])];
        let fake_first_best_bin = Bin {
                bin_contigs: fake_first_best_contigs,
                completeness: 99.9,
                contamination: 0.00,
                bin_type: BinType::prokaryote,
                bin_hash: "fake_hsh_1".to_string()
                };
        let fake_second_best_contigs = vec![Arc::clone(&fake_contigs[2]), Arc::clone(&fake_contigs[4]), Arc::clone(&fake_contigs[5]), Arc::clone(&fake_contigs[6])];
        let fake_second_best_bin = Bin {
                bin_contigs: fake_second_best_contigs,
                completeness: 90.0,
                contamination: 1.00,
                bin_type: BinType::prokaryote,
                bin_hash: "fake_hsh_2".to_string()
                };
        let fake_worst_contigs = vec![Arc::clone(&fake_contigs[0]), Arc::clone(&fake_contigs[1]), Arc::clone(&fake_contigs[5]), Arc::clone(&fake_contigs[6])];
        let fake_worst_bin = Bin {
                bin_contigs: fake_worst_contigs,
                completeness: 1.0,
                contamination: 100.00,
                bin_type: BinType::prokaryote,
                bin_hash: "fake_hsh_3".to_string()
                };
        vec![fake_first_best_bin, fake_second_best_bin, fake_worst_bin]

    } */
    #[test]
    fn run_bin_test() {
        
        let fake_bin_generator = Box::new(FakeBinGen); // fake bin generator works by completeness and contamination based on number of contigs. 
        let fake_contigs = create_fake_contigs_for_unit_test();
        let fake_contig_set_1 = vec![Arc::clone(&fake_contigs[0]), Arc::clone(&fake_contigs[1]), Arc::clone(&fake_contigs[2]), Arc::clone(&fake_contigs[3]), Arc::clone(&fake_contigs[4])]; // this should be the best
        let fake_contig_set_2 = vec![Arc::clone(&fake_contigs[2]), Arc::clone(&fake_contigs[3]), Arc::clone(&fake_contigs[4]), Arc::clone(&fake_contigs[5]), Arc::clone(&fake_contigs[6])]; // this is second best but loses 3 contigs
        let fake_contig_set_3 = vec![Arc::clone(&fake_contigs[6])]; // we should never see set 3 as the only contig is in set 2
        let fake_contig_sets = vec![fake_contig_set_1, fake_contig_set_2, fake_contig_set_3];
        let fake_bin_scorer = BinScorer::initialise_bin_scorer(0.5, 0.5);
        let fake_results = run_classic_best_bin(fake_contig_sets, fake_bin_generator, &fake_bin_scorer);

        assert_eq!(fake_results.len(), 2);

    }
    struct FakeBinGen;

    impl BinGenerator for FakeBinGen {
        fn generate_new_bin_from_contigs(&self, contigs: Vec<Arc<Contig>>) -> Option<Bin> {  // fake bin generator for testing which basically a bins quality is higher if it has more contigs
            let contig_length = contigs.len();
            Some(Bin {
                bin_contigs: contigs,
                completeness: contig_length as f64 / 50.0,
                contamination: contig_length as f64 / 100.0,
                bin_type: BinType::prokaryote,
                bin_hash: format!("fake_{}_sized_bin", contig_length)
            })
    }
}


    
}
