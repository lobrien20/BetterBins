use blake3;
use itertools::Itertools;
use rayon::prelude::{IntoParallelRefIterator, ParallelIterator};
use std::fs::{File, self};
use std::io::{Read, Write};
use std::path::PathBuf;
use std::sync::Arc;
use std::sync::atomic::{AtomicUsize, Ordering};


use crate::contigs::Contig;


pub fn get_contig_headers_from_fasta_file(fasta_path: &PathBuf) -> Vec<String> {
    let mut file = File::open(fasta_path).expect("Error: could not open fasta file");
    let mut contig_file_text = String::new();
    file.read_to_string(&mut contig_file_text);
    let contig_headers: Vec<String> = contig_file_text.lines().filter(|x| x.contains(">")).map(|x| x.to_string()).collect();
    contig_headers
}

pub fn create_new_fasta_file(contig_refs: &Vec<Arc<Contig>>, fasta_path: &PathBuf) {

    let mut fasta_string = String::new();
    let fasta_headers: Vec<String> = contig_refs.iter().map(|x| x.header.clone()).collect();
    let fasta_sequences: Vec<String> = contig_refs.iter().map(|x| x.sequence.clone()).collect();

    for (header, sequence) in fasta_headers.iter().zip(fasta_sequences.iter()) {

        let contig_string = format!("{}\n{}\n", &header, &sequence);
        fasta_string.push_str(&contig_string);

    }

    let mut fasta_file = File::create(fasta_path).expect("can't create fasta file");
    fasta_file.write_all(fasta_string.as_bytes()).expect("Error could not write fasta string to file");

}
pub fn create_bin_fasta(contig_refs: &[Arc<Contig>], fasta_path: &PathBuf) {

    let mut fasta_string = String::new();
    let fasta_headers: Vec<String> = contig_refs.iter().map(|x| x.header.clone()).collect();
    let fasta_sequences: Vec<String> = contig_refs.iter().map(|x| x.sequence.clone()).collect();

    for (header, sequence) in fasta_headers.iter().zip(fasta_sequences.iter()) {

        let contig_string = format!("{}\n{}\n", &header, &sequence);
        fasta_string.push_str(&contig_string);

    }

    let mut fasta_file = File::create(fasta_path).expect(&format!("can't create fasta file at path: {}", &fasta_path.to_str().unwrap()));
    fasta_file.write_all(fasta_string.as_bytes()).expect("Error could not write fasta string to file");

}


pub fn generate_hash_from_contig_comb(contig_combination: &Vec<String>) -> String {

    let contig_combo_stripped: Vec<String> = contig_combination.into_iter().map(|x| x.trim().to_string()).collect();
    let contig_comb_string = contig_combo_stripped.join("");
    let hash = blake3::hash(contig_comb_string.as_bytes());
    let contig_combined_hash = hash.to_hex().to_string();
    contig_combined_hash
}

pub fn generate_hash_from_contigs(contig_combination: &[Arc<Contig>]) -> String {
    
    let contig_combo_string = contig_combination.iter().map(|contig| 
        {contig.header.trim()}).join("");
    
    let hash = blake3::hash(contig_combo_string.as_bytes());
    let contig_combined_hash = hash.to_hex().to_string();
    
    contig_combined_hash
}


pub fn generate_hash_from_contig_arc_comb(contig_combination: &Vec<String>) -> String {

    let contig_combo_stripped: Vec<String> = contig_combination.into_iter().map(|x| x.trim().to_string()).collect();
    let contig_comb_string = contig_combo_stripped.join("");
    let hash = blake3::hash(contig_comb_string.as_bytes());
    let contig_combined_hash = hash.to_hex().to_string();
    contig_combined_hash
}

pub fn check_and_remove_bads_in_hash_directory(hash_directory: &PathBuf) -> AtomicUsize {
    let bin_dirs = fs::read_dir(hash_directory).unwrap();
    let bin_dir_paths: Vec<PathBuf> = bin_dirs.into_iter().map(|x| x.unwrap().path()).collect();
    let unfinished_dirs_count = AtomicUsize::new(0);
    bin_dir_paths.par_iter().for_each(|x| {
            if !x.join("done_file.txt").is_file() {
                if let Err(e) = fs::remove_dir_all(x) {
                    eprintln!("Failed to remove directory: {}", e);
                } else {
                    unfinished_dirs_count.fetch_add(1, Ordering::Relaxed);
                }
            }
        });
    unfinished_dirs_count

}

pub fn check_hash_directory_not_too_big(hash_directory: &PathBuf) {
    let bin_dirs = fs::read_dir(hash_directory).unwrap();
    let bin_dir_paths: Vec<PathBuf> = bin_dirs.into_iter().map(|x| x.unwrap().path()).collect();
    if bin_dir_paths.len() > 1000 {
        panic!("Bin hash directory too big");
    }
}

/* 
pub fn test_permutation_himem(permutation: Vec<Vec<Arc<Contig>>>, bin_generator: Arc<BinGen>, bin_info_storer: Arc<BinInfoStorer>, bin_scorer: Arc<BinScorer>) -> BinSet {

    let mut bin_set_result = BinSet::constructnew();
    let mut used_contigs = Vec::new();

    for contig_set in permutation {
        let available_contigs: Vec<Arc<Contig>> = contig_set.clone().into_iter().filter(|x| !used_contigs.contains(x)).collect();
        if available_contigs.len() == 0 {
            continue
        }

        match bin_info_storer.check_whether_contig_combo_in_hashmap(&available_contigs) {
            Some(res) => {
                if res.1 > bin_generator.minimum_completeness && res.2 < bin_generator.maximum_contamination {
                    used_contigs.extend(available_contigs.clone());
                    bin_set_result.add_new_best_bin((available_contigs, res.0, res.1, res.2));
                }
                continue

            },
            None => { match bin_generator.create_new_bin_from_contigs(&available_contigs) {

            Some(result) => {
                used_contigs.extend(available_contigs.clone());
                bin_set_result.add_new_best_bin((available_contigs, result.0, result.1, result.2))

                },
            None => {}
            }
            }
        }

    }

    bin_set_result
}
*/


mod tests {
    use crate::{contigs::Contig};

    use super::*;
    use std::{fs::{self, DirEntry, remove_dir_all}, env::{self, temp_dir}, path::Path};
    use fs_extra::dir::CopyOptions;
    use lazy_static::lazy_static;
 //   lazy_static! {
   //     static ref TEST_DIR: PathBuf = PathBuf::from("tests/test_data/utils_test/");
   //  }

    #[test]
    fn test_successful_removal_of_unfinished_dirs() {
        let example_dir_path_pathbuf = PathBuf::from("tests/test_data/temp_remove_dir_test/example_directory/");
        fs::create_dir_all("tests/test_data/temp_remove_dir_test/example_directory/example_bad/");
        fs::create_dir_all("tests/test_data/temp_remove_dir_test/example_directory/example_good/");
        File::create("tests/test_data/temp_remove_dir_test/example_directory/example_good/done_file.txt");
        let remaining_number_of_dirs = check_and_remove_bads_in_hash_directory(&example_dir_path_pathbuf);
        assert_eq!(remaining_number_of_dirs.load(Ordering::Relaxed), 1);
    }

}
