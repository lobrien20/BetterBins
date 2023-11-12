use std::{sync::Arc, path::PathBuf, fs::{self, File}, io::Write};

use fs_extra::dir::CopyOptions;
use itertools::Itertools;
use log::debug;

use crate::bin_info_storage::Bin;

#[derive(Clone, Debug, PartialEq)]
pub struct BinSet {
    pub bins: Vec<Arc<Bin>>,
    pub bin_set_order: Option<Vec<usize>>
}

impl BinSet {
    pub fn constructnew() -> BinSet {
        
        BinSet {bins: Vec::new(), bin_set_order: None }
        
    }
    pub fn make_bin_set_from_bins_vec(bins_vec: Vec<Bin>) -> BinSet {
        let bin_arc_vec = bins_vec.into_iter().map(|bin| Arc::new(bin)).collect_vec(); // note, creating arc is only relevant for non classical algorithm, but less boilerplate this way
        BinSet { bins: bin_arc_vec, bin_set_order: None }
    }

    
    pub fn add_bin(&mut self, bin_to_add: Arc<Bin>) {
      
        self.bins.push(bin_to_add);
      
    }

    pub fn create_best_bin_dir_and_info_from_best_hashes(&self, hash_directory: &PathBuf, best_bin_directory: &PathBuf, copy_bins: bool) {
        debug!("Creating directory for best bin at: {:?}", best_bin_directory);

        fs::create_dir(best_bin_directory).unwrap();
        let options = CopyOptions::new();
        let mut best_bin_infos = Vec::new();
        for bin in &self.bins {
            let hash_path = hash_directory.join(&bin.bin_hash);

            if copy_bins == true {

                fs_extra::dir::copy(&hash_path, &best_bin_directory, &options).unwrap();

            }

            best_bin_infos.push((bin.bin_hash.to_string(), bin.bin_type.to_string(), bin.completeness, bin.contamination, bin.bin_contigs.len()));

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




}