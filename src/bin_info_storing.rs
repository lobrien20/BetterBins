use std::{collections::HashMap, path::PathBuf, os, fs::{self, File}, sync::Arc, io::Write, cmp::Ordering};

use blake3::Hash;
use fs_extra::dir::CopyOptions;
use itertools::Itertools;
use log::debug;
use non_dominated_sort::DominanceOrd;
use non_dominated_sort::non_dominated_sort;
use crate::{contigs::Contig, bin_classes::{Bin, Checkm2BinQualityGetter, GetBinQuality, CompleasmBinQualityGetter, BinGenerator}, evo_alg_utils::best_bins_dominance_ord};


pub struct BinScorer {
    pub contamination_weight: f64,
    pub completion_weight: f64,
    pub scoring_approach: BinSetScore,
    pub evo_contamination_scale_power: f64,
    pub evo_completion_scale_power: f64
}

impl BinScorer {

    pub fn score_bin(&self, completion_value: f64, contamination_value: f64) -> f64 {

        let bin_score = ((&self.completion_weight * completion_value) - (&self.contamination_weight * contamination_value)) / 100.0;
        bin_score

    }

    pub fn score_bin_set(&self, bin_set: &BinSet) -> f64 {
        let mut total_score = 0.0;
        for (contigs, string, comp, cont) in bin_set.best_bins.iter() {
            total_score += self.score_bin(*comp, *cont);
        }
        total_score
    }


    pub fn get_best_bin_sets(&self, binsets: &Vec<BinSet>) -> Vec<BinSet> {
        match &self.scoring_approach {
            BinSetScore::CompContWeighted => {
                self.get_best_bin_sets_using_comp_cont_weighted(binsets)
            }
            BinSetScore::SharedCompCont => {
                self.get_best_bin_sets_through_non_dom_sort(binsets)
            }
        }

    }

    fn get_best_bin_sets_using_comp_cont_weighted(&self, binsets:  &Vec<BinSet>) -> Vec<BinSet> {
        let mut ranked_bin_sets = binsets.clone();
        ranked_bin_sets.sort_by(|a, b| b.overall_score.unwrap().partial_cmp(&a.overall_score.unwrap()).unwrap());
        let mut best_sets = vec![ranked_bin_sets[0].clone()];
        ranked_bin_sets.remove(0);
        for binset in ranked_bin_sets {
            if binset.overall_score.unwrap() == best_sets[0].overall_score.unwrap() {
                best_sets.push(binset)
            }
        }
        best_sets
    }

    fn get_best_bin_sets_through_non_dom_sort(&self, binsets: &[BinSet]) -> Vec<BinSet> {
        let mut result = non_dominated_sort(binsets, self);
        let mut last_front = Vec::new();

        while !result.is_empty() {
            last_front = result.iter().map(|x| x.0.clone()).collect();
            result = result.next_front();
        }

        last_front
    }


    pub fn get_least_contam_and_most_complete(&self, binsets: &Vec<BinSet>) -> (BinSet, BinSet){
        if binsets.len() == 1 {
            return (binsets[0].clone(), binsets[0].clone())
        }
        let mut ranked_bin_sets = binsets.clone();
        ranked_bin_sets.sort_by(|a, b| b.overall_completeness.partial_cmp(&a.overall_completeness).unwrap());
        let most_complete = ranked_bin_sets.remove(0);
        ranked_bin_sets.sort_by(|a, b| b.overall_contamination.partial_cmp(&a.overall_contamination).unwrap());
        let least_contam = ranked_bin_sets.pop().unwrap();
        (least_contam, most_complete)
    }



}

pub enum BinSetRankCriteria {
    MostComplete,
    LeastContaminated,
    LeastBins
}




impl DominanceOrd for BinScorer {
    type T = BinSet;

    fn dominance_ord(&self, a: &Self::T, b: &Self::T) -> Ordering {
        best_bins_dominance_ord(a, b)
    }
}


#[derive(Clone, Debug)]
pub enum BinSetScore {
    CompContWeighted,
    SharedCompCont
}


#[derive(Clone, Debug, PartialEq)]
pub struct BinSet {
    pub best_bins: Vec<(Vec<Arc<Contig>>, String, f64, f64)>,
    pub overall_score: Option<f64>,
    pub overall_contamination: f64,
    pub overall_completeness: f64,
    pub bin_set_order: Option<Vec<usize>>,
}

impl BinSet {
    pub fn constructnew() -> BinSet {
        BinSet {
            best_bins: Vec::new(),
            overall_score: None,
            overall_completeness: 0.0,
            overall_contamination: 0.0,
            bin_set_order: None
        }
    }


    pub fn add_new_best_bin(&mut self, bin_to_add: (Vec<Arc<Contig>>, String, f64, f64)) {
  //      self.overall_completeness += bin_to_add.2;

   //     self.overall_contamination += bin_to_add.3;
        self.best_bins.push(bin_to_add);

    }

    pub fn get_unweighted_completion(&self) -> f64 {
        self.best_bins.iter().map(|x| x.2).sum()
    }
    pub fn get_unweighted_contamination(&self) -> f64 {
        self.best_bins.iter().map(|x| x.3).sum()
    }

    pub fn generate_overall_comp_and_cont(&mut self, bin_scorer: &BinScorer) {
        
        self.overall_completeness = self.best_bins.iter().map(|x| x.2.powf(bin_scorer.evo_completion_scale_power)).sum();
        self.overall_contamination = self.best_bins.iter().map(|x| x.3.powf(bin_scorer.evo_completion_scale_power)).sum();

    }

        


    pub fn create_best_bin_dir_and_info_from_best_hashes(&self, hash_directory: &PathBuf, best_bin_directory: &PathBuf, copy_bins: bool) {
        debug!("Creating directory for best bin at: {:?}", best_bin_directory);

        fs::create_dir(best_bin_directory).unwrap();
        let options = CopyOptions::new();
        let mut best_bin_infos = Vec::new();
        for bin in &self.best_bins {
            let hash_path = hash_directory.join(&bin.1);

            if copy_bins == true {

                fs_extra::dir::copy(&hash_path, &best_bin_directory, &options).unwrap();

            }

            let bin_info = self.get_bin_type_and_quality_for_info_file(&hash_path);
            best_bin_infos.push((bin.1.to_string(), bin_info.0, bin_info.1, bin_info.2, bin.0.len()));
        }
        self.write_best_bin_info_file(best_bin_directory, best_bin_infos);

    }

    pub fn write_best_bin_info_file(&self, best_bin_directory: &PathBuf, best_bin_infos: Vec<(String, String, f64, f64, usize)>) {
        let best_bin_file_path = best_bin_directory.join("best_bins_information.tsv");
        let mut best_bin_info_file = File::create(best_bin_file_path).expect("can't create fasta file");
        let best_bin_info_string = best_bin_infos.iter()
        .map(|x| format!("{}\t{}\t{}\t{}\t{}", x.0, x.1, x.2, x.3, x.4)).join("\n");
        let mut best_bin_info_lines = "bin_name\tbin_type\tcompleteness\tcontamination\tnumber_of_contigs\n".to_string();
        best_bin_info_lines.push_str(&best_bin_info_string);

        best_bin_info_file.write_all(best_bin_info_lines.as_bytes()).expect("Error could not write fasta string to file");
    }

    fn get_bin_type_and_quality_for_info_file(&self, bin_dir_path: &PathBuf) -> (String, f64, f64) {
        let prok_qual_report_path = bin_dir_path.join("checkm2_results/quality_report.tsv");
        let euk_qual_report_path =  &bin_dir_path.join("compleasm_results/summary.txt");


        if prok_qual_report_path.is_file() {
            let quality = Checkm2BinQualityGetter::get_comp_and_cont_from_report_file_path(&prok_qual_report_path).unwrap();
            ("prokaryote".to_string(), quality.0, quality.1)

        } else if euk_qual_report_path.is_file() {

            let quality = CompleasmBinQualityGetter::get_comp_and_cont_from_report_file_path(&euk_qual_report_path).unwrap();
            ("eukaryote".to_string(), quality.0, quality.1)

        } else {
            panic!("Couldn't find bin quality report information for best bin at path: {}", &bin_dir_path.to_string_lossy());
        }
    }

}




#[derive(Clone, Debug)]
pub struct BinInfoStorer {

    arc_contig_combo_hashmap: HashMap<Vec<Arc<Contig>>, (String, f64, f64, usize)>,
}

impl BinInfoStorer {
    pub fn initialise_bin_info_storer() -> BinInfoStorer {
        let mut bin_info_storer = BinInfoStorer {
            arc_contig_combo_hashmap: HashMap::new(),
        };
        bin_info_storer
    }



    pub fn add_bins_from_bin_sets_to_arc_contig_hashmap(&mut self, bin_sets: &Vec<BinSet>) {

        bin_sets.iter().flat_map(|bin_set| &bin_set.best_bins).for_each(|x| {
                    self.add_combo_hash_to_arc_hashmap(x.0.clone(), x.1.clone(), (x.2, x.3))});

    }


    pub fn check_whether_contig_combo_in_hashmap(&self, contig_combo: &[Arc<Contig>] ) -> Option<(String, f64, f64, usize)> {

        match self.arc_contig_combo_hashmap.get(contig_combo) {

            Some(quality_score) => return Some(quality_score.clone()),
            None => None

        }
    }


    fn add_combo_hash_to_arc_hashmap(&mut self, contig_combo: Vec<Arc<Contig>>, combo_hash: String, quality_score: (f64, f64)) {
        let number_of_contigs = contig_combo.len();
        self.arc_contig_combo_hashmap.insert(contig_combo, (combo_hash, quality_score.0, quality_score.1, number_of_contigs));

        }

}




mod tests {

    use super::*;
    use std::{fs, env, path::Path};
    use lazy_static::lazy_static;
    lazy_static! {
    static ref TEST_DATA_HASH: PathBuf = PathBuf::from("tests/test_data/bin_classes_testing/hash_results/");
}
    lazy_static! {
    static ref TEST_DATA_DIR: PathBuf = PathBuf::from("tests/test_data/full_test/temp_full_test_run/hash_results/");
}

    #[test]
    fn test_search_for_bins_in_hash_dir() { // assumes that hash results from integration testing is there (it should be)
        let bin_info_storer = BinInfoStorer::initialise_bin_info_storer();
    }

    fn create_some_test_permutations() -> Vec<Vec<usize>> {

        let permut_1 = vec![0, 1, 2, 3, 4, 6, 5];
        let permut_2 = vec![2, 4, 3, 1, 5, 0, 6];
        let permut_3 = vec![6, 4, 1, 5, 3, 0, 2];
        let permut_4 = vec![1, 4, 2, 5, 0, 3, 6];

        vec![permut_1, permut_2, permut_3, permut_4]
    }

    fn create_some_test_best_binset() -> Vec<BinSet> {
        // creates four test best bins.
        // overall score has a difference of 10 between each one (for weighted picking)
        // overall completeness has the same, but the third best bin has a higher contamination than the rest.
        let test_permutations = create_some_test_permutations();
        let mut test_best_bins = Vec::new();
        let mut score_to_add = 0.0;

        for permut in test_permutations {

            let mut best_bins_t= BinSet::constructnew();

            best_bins_t.overall_completeness = 50.0;
            best_bins_t.overall_contamination = 0.0;
            best_bins_t.overall_score = Some(score_to_add);
            best_bins_t.bin_set_order = Some(permut);
            test_best_bins.push(best_bins_t);
            score_to_add += 10.0;


        }
        test_best_bins[3].overall_contamination = 10.0;
        test_best_bins[2].overall_completeness = 90.0;
        test_best_bins
    }
    #[test]
    fn test_score_bin_set() {
        let test_binset = create_some_test_best_binset();
        let bin_scorer = BinScorer {completion_weight: 0.5, contamination_weight: 0.5, scoring_approach:BinSetScore::CompContWeighted, evo_completion_scale_power: 1.0, evo_contamination_scale_power: 1.0};
        let fake_bin_set = vec![(Vec::new(), "na".to_string(), 50.0, 10.0), (Vec::new(), "na".to_string(), 90.0, 5.0)];
        let fake_BinSet = BinSet {best_bins: fake_bin_set, overall_score: None, overall_completeness: 0.0, overall_contamination: 0.0, bin_set_order: None};
        let result = fake_BinSet.get_unweighted_completion();
        assert_eq!(140.0, result);
        let result2 = bin_scorer.score_bin_set(&fake_BinSet);
        assert_eq!(0.625, result2);
        let best_result = bin_scorer.get_best_bin_sets_through_non_dom_sort(&test_binset);
        println!("{:?}", best_result);
    }

}
