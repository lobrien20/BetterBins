use std::{path::PathBuf, fs::{File, self}, io::{Read, Write}, sync::Arc, num, collections::{HashSet, HashMap}, thread, time::Duration, hash};

use itertools::Itertools;
use log::info;
use rayon::prelude::{IntoParallelRefIterator, ParallelIterator, IntoParallelIterator};

use crate::{contigs::{Contig, ContigType, EukaryoticContigInformation}, prokaryotic_contig_gatherer::ProkaryoticBinQualityGetter, eukaryotic_contig_gatherer::EukaryoticBinQualityGetter, contig_type_predictor::ContigTypePredictor, utils::{generate_hash_from_contigs, create_new_fasta_file, create_bin_fasta}, bin_info_storage::{BinType, BinInfoStorage, Bin}, bin_generator::{BinGen, BinTypePrediction, BinGenerator}};

pub fn initialise_tool_through_getting_original_bins_and_contigs(contig_file: &PathBuf, output_directory: &PathBuf, prokaryote_db_path: String, 
    threads: usize, bin_directory_path: &PathBuf, eukaryota_db_path: String, number_of_markers: usize, compleasm_db_type: String, hash_directory_path: &PathBuf, 
    maximum_contamination: f64, minimum_completeness: f64, bin_type_predictor: Box<dyn BinTypePrediction>, bin_info_storer: BinInfoStorage) -> (BinGen, Vec<Bin>) {

    let prok_bin_getter = ProkaryoticBinQualityGetter::initialise(&prokaryote_db_path);
    let euk_bin_getter = EukaryoticBinQualityGetter::initialise(&eukaryota_db_path, number_of_markers, compleasm_db_type);
    let all_contigs = gather_contigs_and_add_info(&contig_file, &output_directory.join("contig_outputs/"), &prok_bin_getter, &euk_bin_getter, threads);
    let bin_generator = BinGen::initialise(Some(prok_bin_getter), Some(euk_bin_getter), hash_directory_path.clone(), maximum_contamination, minimum_completeness, 
    bin_info_storer, bin_type_predictor);
    
    let fasta_files = find_bin_fastas(&bin_directory_path);
    let path_to_bin_hashmap = analyse_initial_bin_fastas_and_generate_path_to_bin_dict(fasta_files, all_contigs, &bin_generator);
    create_summary_file_of_initial_bin_fastas(&path_to_bin_hashmap, output_directory, &bin_directory_path);
    let the_bins = path_to_bin_hashmap.into_values().filter_map(|bin| bin).collect_vec();
    (bin_generator, the_bins)
}




pub fn gather_contigs_and_add_info(contig_file: &PathBuf, output_directory: &PathBuf, 
    prok_bin_getter: &ProkaryoticBinQualityGetter, euk_bin_getter: &EukaryoticBinQualityGetter, threads: usize) -> Vec<Arc<Contig>> {

    let mut all_contigs: Vec<Contig> = generate_all_contigs_from_contig_file(&contig_file);
    fs::create_dir(output_directory);

    ContigTypePredictor::predict_contig_types(&contig_file, &mut all_contigs, output_directory).unwrap();
    prok_bin_getter.add_prok_info_to_contigs_using_checkm2(&contig_file, output_directory, &mut all_contigs, threads);
    euk_bin_getter.add_euk_info_to_contigs_using_compleasm(&contig_file, &mut all_contigs, output_directory, threads);
    
    let all_arc_contigs: Vec<Arc<Contig>> = all_contigs.into_iter().map(|contig| Arc::new(contig)).collect();
    all_arc_contigs
}



pub fn generate_all_contigs_from_contig_file(contig_file: &PathBuf) -> Vec<Contig> {   
      
    let fasta_info = get_fasta_info_from_file(&contig_file);
    let mut bin_contigs = Vec::new();
    
    for header_and_seq in fasta_info {
        
        let contig_fasta_lines: Vec<&str> = header_and_seq.lines().collect();
        let contig_header = contig_fasta_lines[0].trim().to_string();
        let contig_fasta_sequence: String = contig_fasta_lines.into_iter().skip(1).collect::<Vec<_>>().join("");
        let contig = Contig::new_contig(contig_header, contig_fasta_sequence);
        bin_contigs.push(contig);
    }


    bin_contigs
    

}


fn get_fasta_info_from_file(fasta_file_path: &PathBuf) -> Vec<String> {
        
    let mut fasta_file = File::open(&fasta_file_path).expect(&format!("Could not open fasta file for bin at path: {}", &fasta_file_path.to_string_lossy()));
    let mut fasta_lines = String::new();
    fasta_file.read_to_string(&mut fasta_lines).expect(&format!("Could not read fasta file to string"));
    let mut fasta_vec: Vec<String> = fasta_lines.split(">").map(|x| x.to_string()).collect();
    let fasta_info: Vec<String> = fasta_vec.iter().skip(1).map(|x| format!(">{}", x)).collect();
    fasta_info

}



fn find_bin_fastas(bin_dir: &PathBuf) -> Vec<PathBuf> {
    let mut fasta_paths = Vec::new();
    let mut found_bin_in_dir = false;
    for entry in fs::read_dir(bin_dir).unwrap() {
        let entry = entry.unwrap();
        let data = entry.metadata().unwrap();
        let path = entry.path();
        if data.is_file() {
            if let Some(ext) = path.extension() {
                if ext == "fa" || ext == "fasta" {
                    found_bin_in_dir = true;
                    fasta_paths.push(path)
                }
            }
        } else if data.is_dir() {
            // recursively looks for other fastas
            fasta_paths.extend(find_bin_fastas(&path));
        }

    }
    if found_bin_in_dir == true {
        info!("Found bins in directory path: {}", &bin_dir.to_string_lossy());
    }
    fasta_paths
}

fn analyse_initial_bin_fastas_and_generate_path_to_bin_dict(bin_fastas: Vec<PathBuf>, all_contigs: Vec<Arc<Contig>>, bin_gen: &BinGen) -> HashMap<PathBuf, Option<Bin>>{
    println!("found {} fastas", bin_fastas.len());
    let initial_path_to_generated_bin_dict: HashMap<PathBuf, Option<Bin>> = bin_fastas.into_par_iter().map(|fasta_path| {

        let bin_contigs = get_contig_set_arc_from_bin_fasta_and_contigs(&fasta_path, &all_contigs);
        (fasta_path, bin_gen.generate_new_bin_from_contigs(bin_contigs))

    }).collect();
    initial_path_to_generated_bin_dict
}

fn get_contig_set_arc_from_bin_fasta_and_contigs(fasta_file_path: &PathBuf, all_contigs: &Vec<Arc<Contig>>) -> Vec<Arc<Contig>> {
    
    let mut fasta_file = File::open(&fasta_file_path).expect(&format!("Could not open fasta file for bin at path: {}", &fasta_file_path.to_string_lossy()));
    let mut fasta_file_string = String::new();
    fasta_file.read_to_string(&mut fasta_file_string).expect(&format!("Could not read fasta file to string"));
    let fasta_lines = fasta_file_string.lines();
    
    let headers: Vec<&str> = fasta_lines.filter(|line| line.contains(">")).map(|header| header.trim()).collect();
    let mut bin_contigs = Vec::new();
    for header in headers {
        let mut contig_found = false;

        for contig in all_contigs {
            
            if contig.header == header {
            
                bin_contigs.push(Arc::clone(&contig));
                contig_found = true;
                break
            
            }
        }
        if contig_found == false {
            panic!("Error, found a contig in a bin not found in sample contigs file")
        }

    }
    println!("Bin has {} contigs", bin_contigs.len());
    bin_contigs

}


fn create_summary_file_of_initial_bin_fastas(path_to_bin_hashmap: &HashMap<PathBuf, Option<Bin>>, output_directory: &PathBuf, bin_directory: &PathBuf) {
    
    let summary_file_path = output_directory.join("original_bins_summary_file.tsv");
    let mut summary_file_string = "bin_name\tsub_directory\tbin_type\tcompleteness\tcontamination\tnumber_of_contigs\n".to_string();
    let mut summary_file = File::create(summary_file_path).expect("can't create fasta file");
    let mut bin_directory_name = bin_directory.file_name().unwrap();
    
    for (path, bin) in path_to_bin_hashmap {
        let mut directory_of_individual_bin = path.parent().unwrap();
        let mut bin_sub_directory = path.parent().unwrap().file_name().unwrap().to_str().unwrap();
        
        if directory_of_individual_bin == bin_directory_name {
            bin_sub_directory = "na";
        } 
        let bin_name = path.file_name().unwrap().to_str().unwrap();
        let bin_info_string = match bin {
            Some(bin) => create_bin_info_string(bin_name, bin_sub_directory, &bin),
            None => format!("{}\t{}\tna\tna\tna\tna\n", &bin_name, bin_sub_directory)
        };
        summary_file_string.push_str(&bin_info_string);

    }
    info!("The following bins were generated: {}", &summary_file_string);
    summary_file.write_all(summary_file_string.as_bytes()).expect("Error could not write fasta string to file");

} 

fn create_bin_info_string(bin_name: &str, bin_sub_directory: &str, bin: &Bin) -> String {
    let bin_type_string = match bin.bin_type {
        BinType::eukaryote => "eukaryote".to_string(),
        BinType::prokaryote => "prokaryote".to_string(),
    };

    let bin_info_string = format!("{}\t{}\t{}\t{}\t{}\t{}\n", &bin_name, bin_sub_directory, bin_type_string, bin.completeness, bin.contamination, bin.bin_contigs.len());
    bin_info_string

}




#[cfg(test)]
mod tests {
    use crate::{contigs::{ProkaryoticContigInformation, ContigType}, initialise_bins_and_contigs::generate_all_contigs_from_contig_file, bin_generator::MinimumEukMarkerGenes};

    use super::*;
    use std::{fs, env, path::Path};
    use itertools::Itertools;
    use lazy_static::lazy_static;
    use regex::Regex;
    lazy_static! {
    static ref TEST_DATA_DIR: PathBuf = PathBuf::from("tests/new_tests/test_bins/");
    static ref TEST_CONTIGS_PATH: PathBuf = PathBuf::from("tests/new_tests/all_unique_contigs.fa");
    static ref INITIALISE_BC_PATH: PathBuf = PathBuf::from("tests/new_tests/initialise_bins_and_contigs_test/");
}

    lazy_static! { // unit test will only work if database is added here, not uploaded to git repo due to size
        static ref CHECKM2_DB_PATH: PathBuf = PathBuf::from("tests/new_tests/databases_for_testing/uniref100.KO.1.dmnd");
        static ref COMPLEASM_DB_LIB: PathBuf = PathBuf::from("tests/new_tests/databases_for_testing/eukaryota_odb10");

    }
    #[test]
    fn test_initialise_tool() {
        let bin_type_predictor = MinimumEukMarkerGenes { minimum_marker_gene_count: 128 };  
        let contig_file_path = &TEST_DATA_DIR.join("all_unique_contigs.fa");
        let full_initialise_test = &INITIALISE_BC_PATH.join("full_initialise_test");
        let test_bin_dir = &INITIALISE_BC_PATH.join("test_bins/");
        let hash_dir_path = &INITIALISE_BC_PATH.join("hash_test_dir/");
        let bin_info_storer = BinInfoStorage::initialise_bin_info_storer();
        let (bin_gen, bins) = initialise_tool_through_getting_original_bins_and_contigs(contig_file_path, full_initialise_test, CHECKM2_DB_PATH.clone().into_os_string().into_string().unwrap(), 6, &TEST_DATA_DIR, 
            COMPLEASM_DB_LIB.clone().into_os_string().into_string().unwrap(), 
            255, "eukaryota_odb10".to_string(), hash_dir_path, 
            100.0, 0.0, Box::new(bin_type_predictor), bin_info_storer);
            assert_eq!(bins.len(), 7);
            let mut all_bin_completenesses = bins.iter().map(|bin| bin.completeness).collect_vec();
            let mut all_bin_contaminations = bins.iter().map(|bin| bin.contamination).collect_vec();
            all_bin_completenesses.sort_by(|a, b| a.partial_cmp(b).unwrap());
            all_bin_contaminations.sort_by(|a, b| a.partial_cmp(b).unwrap());
            let expected_bin_contaminations = vec![0.0, 0.0, 8.78, 8.78, 8.84, 13.28, 13.28];
            let expected_bin_completnesses = vec![52.9, 52.9, 92.94, 93.33, 100.0, 100.0, 100.0];
            assert_eq!(expected_bin_contaminations, all_bin_contaminations);
            assert_eq!(expected_bin_completnesses, all_bin_completenesses);
    }



    #[test]
    fn test_gather_contigs_and_get_info() {
        let contig_file_path = &TEST_DATA_DIR.join("all_unique_contigs.fa");
        let prok_bin_getter = ProkaryoticBinQualityGetter::initialise(&CHECKM2_DB_PATH.to_string_lossy());
        let euk_bin_getter = EukaryoticBinQualityGetter::initialise(&COMPLEASM_DB_LIB.to_string_lossy(), 255, "eukaryota_odb10".to_string());
        let contigs = gather_contigs_and_add_info(contig_file_path, &INITIALISE_BC_PATH, &prok_bin_getter, &euk_bin_getter, 6);
        assert_eq!(contigs.len(), 1504);
        let expected_diamond_info_path = &INITIALISE_BC_PATH.join("checkm2_results/diamond_output/DIAMOND_RESULTS.tsv");
        let expected_protein_info_path = &INITIALISE_BC_PATH.join("checkm2_results/protein_files/all_unique_contigs.faa");
        let expected_eukrep_file = &INITIALISE_BC_PATH.join("all_contig_eukrep_results.fa");
        let expected_prok_file = &INITIALISE_BC_PATH.join("all_contig_prok_eukrep_results.fa");
        let expected_busco_path = &INITIALISE_BC_PATH.join("compleasm_results/eukaryota_odb10/full_table_busco_format.tsv");
        let expected_files = vec![expected_eukrep_file, expected_prok_file, expected_diamond_info_path, expected_protein_info_path, expected_busco_path];
        for file in expected_files {
            
            if !file.is_file() {
                panic!("test_gather_contigs_and_get_info: {} not produced", &file.to_string_lossy());
            }
            
        }
    }

    fn get_unique_diamonds_for_test(prok_info: &ProkaryoticContigInformation) -> Vec<&str>{
        let diamond_vec = prok_info.diamond_lines.iter()
                    .map(|line| {
                        let re = Regex::new(r"(NODE_[^_]+_length_[^_]+_cov_[^_]+)").unwrap();
                        let caps = re.captures(line).unwrap();
                    
                        // Extract the captured group which excludes the trailing _number
                        caps.get(0).map_or("", |m| m.as_str())
                    }).unique().collect_vec();
      // println!("{:?}", diamond_vec);
        diamond_vec
    }
}