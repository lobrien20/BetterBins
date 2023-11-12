use std::{path::PathBuf, io::{self, Read}, fs::{self, File}, sync::Arc, collections::HashMap};
use blake3::Hash;
use itertools::Itertools;
use log::{info, debug};
use rayon::{iter::IntoParallelRefIterator, prelude::ParallelIterator};
use crate::{bin_classes::{BinGenerator, Bin, BinType}, contigs::Contig, utils::{get_contig_headers_from_fasta_file, generate_hash_from_contig_comb}};
use glob::glob;

pub fn gather_initial_bins_and_contig_information(bin_dir: &PathBuf, bin_generator_arc: Arc<BinGenerator>) -> (Vec<Bin>, Vec<Arc<Contig>>) {
    
    let mut constructed_bins: Vec<Bin> = initial_bin_search_and_quality_scoring(bin_dir, bin_generator_arc);
    debug!("Found {} bins", constructed_bins.len());
    let bins_mut: Vec<&mut Bin> = constructed_bins.iter_mut().collect();
    info!("Found {} bins", bins_mut.len());
    let unique_contigs = BinContigInterface::generate_unique_contigs_and_add_contig_ref_to_bins(bins_mut);
    info!("Found {} unique contigs", unique_contigs.len());

    (constructed_bins, unique_contigs)
}

pub fn initial_bin_search_and_quality_scoring(bin_dir: &PathBuf, bin_generator_arc: Arc<BinGenerator>) -> Vec<Bin> {


    let mut fasta_paths = find_bin_fastas(bin_dir);
    fasta_paths = ensure_no_duplicate_fastas(fasta_paths);

    let bin_results: Vec<Option<Bin>> = fasta_paths.par_iter().map(|fasta_path| {
        let bin_gen_clone = Arc::clone(&bin_generator_arc);
        println!("creating bin!");
        bin_gen_clone.create_bin_from_original_bin_fasta_file(&fasta_path)
        
        }
    ).collect();
    let constructed_bins: Vec<Bin> = bin_results.into_iter().filter_map(|x| x).collect();

    constructed_bins
}

fn ensure_no_duplicate_fastas(fasta_paths: Vec<PathBuf>) -> Vec<PathBuf> {
    let mut unique_bins = HashMap::new();
    for fasta_path in fasta_paths {
        let contig_headers = get_contig_headers_from_fasta_file(&fasta_path);
        let contig_hash = generate_hash_from_contig_comb(&contig_headers);
        match unique_bins.get(&contig_hash) {
            Some(_) => info!("Found duplicate bin!"),
            None => {
                unique_bins.insert(contig_hash, fasta_path);
            }
        }
    }
    let mut unique_fastas: Vec<PathBuf> = unique_bins.into_values().map(|x| x.clone()).collect();
    unique_fastas
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

pub struct BinContigInterface;

// -> HashMap<Contig, Vec<&'static Bin>>
impl BinContigInterface {
    pub fn generate_unique_contigs_and_add_contig_ref_to_bins(bins: Vec<&mut Bin>) -> Vec<Arc<Contig>> {
        let mut used_contig_hashmap: HashMap<String, Arc<Contig>> = HashMap::new();
        let mut used_contigs = Vec::new();
        for bin in bins {
            let bin_contigs = BinContigInterface::generate_contigs_from_bin(&bin);

            for contig in &bin_contigs {
           
                if let Some(arc_contig_ref) = used_contig_hashmap.get_mut(&contig.header) {
            
                    bin.add_contig_ref_to_bin(Arc::clone(&arc_contig_ref));

                } else {
                   
                    let the_contig = contig.clone();
                    let arc_contig_ref = Arc::new(the_contig);
                    bin.add_contig_ref_to_bin(Arc::clone(&arc_contig_ref));
                    used_contig_hashmap.insert(contig.header.clone(), arc_contig_ref.clone());
                    used_contigs.push(Arc::clone(&arc_contig_ref));
                    
                }
            }

        }

        used_contigs
    }

    pub fn generate_contigs_from_bin(bin: &Bin) -> Vec<Contig> {   
      
        let fasta_info = BinContigInterface::get_fasta_info_from_file(&bin.fasta_path);
        let mut bin_contigs = Vec::new();
        
        for header_and_seq in fasta_info {
            
            let contig_fasta_lines: Vec<&str> = header_and_seq.lines().collect();
            let contig_header = contig_fasta_lines[0].trim().to_string();
            let contig_fasta_sequence: String = contig_fasta_lines.into_iter().skip(1).collect::<Vec<_>>().join("");
            let contig = Contig::new_contig(contig_header, contig_fasta_sequence);
            bin_contigs.push(contig);
        }

        match &bin.bin_type {
            Some(bin_type) => {
                if bin_type == &BinType::prokaryote {
                    let contigs_mut: Vec<&mut Contig> = bin_contigs.iter_mut().collect();
                    BinContigInterface::add_orf_protein_contig_information_for_prok_bin(bin, contigs_mut);
                }
            },
            None => {
                panic!("Error, no bin type for bin when contig generated")
            }
        }

        bin_contigs
        
 
    }
    fn add_orf_protein_contig_information_for_prok_bin(prok_bin: &Bin, contigs_in_bin: Vec<&mut Contig> ) {
        
        let protein_info = BinContigInterface::get_protein_info_from_file(&PathBuf::from(&format!("{}/checkm2_results/protein_files/{}.faa", &prok_bin.bin_dir_path.to_string_lossy(), &prok_bin.bin_hash)));
        let diamond_info = BinContigInterface::get_prokaryotic_diamond_info_from_file(&PathBuf::from(&format!("{}/checkm2_results/diamond_output/DIAMOND_RESULTS.tsv", &prok_bin.bin_dir_path.to_string_lossy())));     
       
        for contig in contigs_in_bin {

            let mut header_without_shark_symbol = contig.header.clone();
            header_without_shark_symbol.remove(0);
            let contig_protein_info: Vec<String> = protein_info.iter().filter(|x| x.contains(&contig.header)).map(|x| x.to_string()).collect();
            let diamond_protein_info: Vec<String> = diamond_info.iter().filter(|x| x.contains(&header_without_shark_symbol)).map(|x| x.to_string()).collect();
            contig.add_diamond_and_protein_information(diamond_protein_info, contig_protein_info);
        }


    }


    
    fn get_fasta_info_from_file(fasta_file_path: &PathBuf) -> Vec<String> {
        
        let mut fasta_file = File::open(&fasta_file_path).expect(&format!("Could not open fasta file for bin at path: {}", &fasta_file_path.to_string_lossy()));
        let mut fasta_lines = String::new();
        fasta_file.read_to_string(&mut fasta_lines).expect(&format!("Could not read fasta file to string"));
        let mut fasta_vec: Vec<String> = fasta_lines.split(">").map(|x| x.to_string()).collect();
        let fasta_info: Vec<String> = fasta_vec.iter().skip(1).map(|x| format!(">{}", x)).collect();
        fasta_info

    }

    fn get_protein_info_from_file(protein_file_path: &PathBuf) -> Vec<String> {
    
        let mut protein_file = File::open(&protein_file_path).expect(&format!("Could not open protein file for bin at path: {}", &protein_file_path.to_string_lossy()));
        let mut protein_lines = String::new();
        protein_file.read_to_string(&mut protein_lines).expect("Could not read protein file to string");
        let mut protein_vec: Vec<String> = protein_lines.split(">").map(|x| x.to_string()).collect();
//        println!("1 is: {}", protein_vec[0]);
   //     println!("2 is: {}", protein_vec[1]);
        let protein_info: Vec<String> = protein_vec.iter().skip(1).map(|x| format!(">{}", x)).collect();
        protein_info
    }


    fn get_prokaryotic_diamond_info_from_file(diamond_file_path: &PathBuf) -> Vec<String> {
    
        let mut diamond_file = File::open(&diamond_file_path).expect(&format!("Could not open diamond file for bin at path: {}", &diamond_file_path.to_string_lossy()));
        let mut diamond_lines = String::new();
        diamond_file.read_to_string(&mut diamond_lines).expect("Could not read diamond file to string");
        let diamond_info: Vec<String> = diamond_lines.lines()
        .map(|x| x.to_string().split("Ω").nth(1).unwrap().to_string())
        .collect();
        diamond_info   

    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use std::{fs, env, path::Path};
    use lazy_static::lazy_static;
    lazy_static! {
    static ref TEST_DATA_HASH: PathBuf = PathBuf::from("tests/test_data/initial_bins_and_contigs_test_data/hash_results/");
    }
    lazy_static! {
        static ref TEST_DATA_DIR: PathBuf = PathBuf::from("tests/test_data/initial_bins_and_contigs_test_data/");
    }
    lazy_static! { // unit test will only work if database is added here, not uploaded to git repo due to size
        static ref COMPLEASM_DB_LIB: PathBuf = PathBuf::from("tests/test_data/databases_for_testing/");
    }
    lazy_static! { // unit test will only work if database is added here, not uploaded to git repo due to size
    
        static ref CHECKM2_DB_PATH: PathBuf = PathBuf::from("tests/test_data/databases_for_testing/uniref100.KO.1.dmnd");
    
    }


    fn create_bin_generator_for_testing() -> BinGenerator {
        BinGenerator::initialise_bin_gen(
                Some(CHECKM2_DB_PATH.to_path_buf()), Some(COMPLEASM_DB_LIB.to_path_buf()), TEST_DATA_HASH.to_path_buf(), 0.0, 0.0)
                
        
    }

    fn create_premade_bin_1() -> Bin {
        // 1284 contigs
        let premade_bin = Bin {
            fasta_path: TEST_DATA_HASH.join("08c929dd57d754fe8e56f2befc8db4daed8bded44d54525e970e9be228345b9c/08c929dd57d754fe8e56f2befc8db4daed8bded44d54525e970e9be228345b9c.fa"),
            bin_contigs: None,
            completeness: Some(51.74),
            contamination: Some(13.2),
            bin_type: Some(BinType::prokaryote),
            bin_hash: "08c929dd57d754fe8e56f2befc8db4daed8bded44d54525e970e9be228345b9c".to_string(),
            bin_dir_path: TEST_DATA_HASH.join("08c929dd57d754fe8e56f2befc8db4daed8bded44d54525e970e9be228345b9c/")
        };
        premade_bin
    }

    fn create_premade_bin_2() -> Bin {
        // 203 contigs
        let mut premade_bin_2 = Bin {fasta_path: TEST_DATA_HASH.join("20ade716c5c524f277bb6a0cccc80e289494a5060df506f1090cf3736ac38b02/20ade716c5c524f277bb6a0cccc80e289494a5060df506f1090cf3736ac38b02.fa"),
        bin_contigs: None,
        completeness: Some(100.0),
        contamination: Some(8.75),
        bin_type: Some(BinType::prokaryote),
        bin_hash: "20ade716c5c524f277bb6a0cccc80e289494a5060df506f1090cf3736ac38b02".to_string(),
        bin_dir_path: TEST_DATA_HASH.join("20ade716c5c524f277bb6a0cccc80e289494a5060df506f1090cf3736ac38b02/")

        };
        premade_bin_2

    }

    
    #[test]
    fn verify_fasta_finding() {
        let bin_fastas = find_bin_fastas(&TEST_DATA_DIR.join("bins/"));
        assert_eq!(bin_fastas.len(), 4);
    }
    #[test]
    fn verify_initial_bin_search() {
       let bin_gener = Arc::new(create_bin_generator_for_testing());

        let bins = initial_bin_search_and_quality_scoring(&TEST_DATA_DIR, bin_gener);
        assert_eq!(bins.len(), 4)
    }
    fn verify_correct_prokaryote_specific_info_obtained() {
        
        let mut premade_bin_to_be_created = Bin {fasta_path: TEST_DATA_HASH.join("20ade716c5c524f277bb6a0cccc80e289494a5060df506f1090cf3736ac38b02/20ade716c5c524f277bb6a0cccc80e289494a5060df506f1090cf3736ac38b02.fa"),
        bin_contigs: None,
        completeness: Some(100.0),
        contamination: Some(8.75),
        bin_type: Some(BinType::prokaryote),
        bin_hash: "20ade716c5c524f277bb6a0cccc80e289494a5060df506f1090cf3736ac38b02".to_string(),
        bin_dir_path: TEST_DATA_HASH.join("/20ade716c5c524f277bb6a0cccc80e289494a5060df506f1090cf3736ac38b02/")

        };

    }
    #[test]
    fn verify_getting_protein_info_from_file() {
        let protein_info_pathbuf = TEST_DATA_HASH.join("20ade716c5c524f277bb6a0cccc80e289494a5060df506f1090cf3736ac38b02/checkm2_results/protein_files/20ade716c5c524f277bb6a0cccc80e289494a5060df506f1090cf3736ac38b02.faa");

        let protein_info = BinContigInterface::get_protein_info_from_file(&protein_info_pathbuf);
        assert_eq!(protein_info.len(), 4872);

        for info in protein_info {
    
            if info.matches(">").count() != 1 {
                println!("{}", info); 
                panic!("Bad shark symbol count in verifying getting protein info from file");

       
            }

        }

    }
    #[test]
    fn verify_getting_diamond_info_from_file() {
        let diamond_info_pathbuf = TEST_DATA_HASH.join("20ade716c5c524f277bb6a0cccc80e289494a5060df506f1090cf3736ac38b02/checkm2_results/diamond_output/DIAMOND_RESULTS.tsv");
        let diamond_info = BinContigInterface::get_prokaryotic_diamond_info_from_file(&diamond_info_pathbuf);
        assert_eq!(diamond_info.len(), 2365);
        for info in diamond_info {

            let diamond_sections = info.split("\t").collect::<Vec<&str>>();
            for sect in &diamond_sections {
                println!("{}", sect);
            }
            if diamond_sections.len() != 12 {
                println!("{}", info);
                panic!("Bad number of sections in a line for verify getting diamond info from file");
            }

        }
    }
    #[test]
    fn verify_prokaryotic_contigs_correctly_constructed() {
        let mut premade_bin_to_be_created = Bin {fasta_path: TEST_DATA_HASH.join("20ade716c5c524f277bb6a0cccc80e289494a5060df506f1090cf3736ac38b02/20ade716c5c524f277bb6a0cccc80e289494a5060df506f1090cf3736ac38b02.fa"),
        bin_contigs: None,
        completeness: Some(100.0),
        contamination: Some(8.75),
        bin_type: Some(BinType::prokaryote),
        bin_hash: "20ade716c5c524f277bb6a0cccc80e289494a5060df506f1090cf3736ac38b02".to_string(),
        bin_dir_path: TEST_DATA_HASH.join("20ade716c5c524f277bb6a0cccc80e289494a5060df506f1090cf3736ac38b02/")

        };
        let mut the_contigs = BinContigInterface::generate_contigs_from_bin(&premade_bin_to_be_created);
        let contigs_mut: Vec<&mut Contig> = the_contigs.iter_mut().collect();
        BinContigInterface::add_orf_protein_contig_information_for_prok_bin(&premade_bin_to_be_created, contigs_mut);
        for contig in the_contigs {
            let prok_contig = contig.prokaryotic_contig_info.unwrap();

            verify_diamond_info_correctly_constructed(&prok_contig.diamond_lines);
            verify_protein_info_correctly_constructed(&prok_contig.predicted_protein_lines);

        }
    }

    fn verify_diamond_info_correctly_constructed(diamond_info: &Vec<String>) {
        for info in diamond_info {

            let diamond_sections = info.split("\t").collect::<Vec<&str>>();
            for sect in &diamond_sections {
                println!("{}", sect);
            }
            if diamond_sections.len() != 12 {
                println!("{}", info);
                panic!("Bad number of sections in a line for verify getting diamond info from file");
            }
        }
    }

    fn verify_protein_info_correctly_constructed(protein_info: &Vec<String>) {
        for info in protein_info {
    
            if info.matches(">").count() != 1 {
                println!("{}", info); 
                panic!("Bad shark symbol count in verifying getting protein info from file");

       
            }

        }
    }
    #[test]
    
    fn verify_getting_fasta_info_from_file() {
        let fasta_info_pathbuf = TEST_DATA_HASH.join("20ade716c5c524f277bb6a0cccc80e289494a5060df506f1090cf3736ac38b02/20ade716c5c524f277bb6a0cccc80e289494a5060df506f1090cf3736ac38b02.fa");
        let fasta_info = BinContigInterface::get_fasta_info_from_file(&fasta_info_pathbuf);
        assert_eq!(fasta_info.len(), 203);

        for info in fasta_info {
    
            if info.matches(">").count() != 1 {
                println!("{}", info); 
                panic!("Bad shark symbol count in verifying getting fasta info from file");

       
            }

        }
    }
    #[test]
    fn validate_construction_of_vector_of_contigs() {
        
        let premade_bin_1 = create_premade_bin_1();
        let bin_1_contigs = BinContigInterface::generate_contigs_from_bin(&premade_bin_1);
        assert_eq!(bin_1_contigs.len(), 1284);
        let premade_bin_2 = create_premade_bin_2();
        let bin_2_contigs = BinContigInterface::generate_contigs_from_bin(&premade_bin_2);
        assert_eq!(bin_2_contigs.len(), 203);
        
    }
    #[test]

    fn ensure_contigs_correctly_added_to_bin() {
        let mut premade_bin_1 = create_premade_bin_1();
        let mut premade_bin_2 = create_premade_bin_2();
        let bins_vec = vec![&mut premade_bin_1, &mut premade_bin_2];
        let unique_contigs = BinContigInterface::generate_unique_contigs_and_add_contig_ref_to_bins(bins_vec);
        assert_eq!(1487, unique_contigs.len());
        match premade_bin_1.bin_contigs {
            Some(contigs) => 
            {
                assert_eq!(contigs.len(), 1284);
                for contig in contigs {
                    if contig.prokaryotic_contig_info == None {
                        panic!("prokaryote contigs not successfully had predicted protein lines added");
                    }
                }
            
            },
            None => panic!("No contigs")
        }
        match premade_bin_2.bin_contigs {
            Some(contigs) => assert_eq!(contigs.len(), 203),
            None => panic!("No contigs")
        }
    }
    #[test]
    fn test_whole_module_through_gather_initial_bins_and_contig_information_method_with_premade_bins() {
        let bin_gener = Arc::new(create_bin_generator_for_testing());
        let bin_gen_ref = Arc::clone(&bin_gener);
        let bins = initial_bin_search_and_quality_scoring(&TEST_DATA_DIR, bin_gen_ref);
        let (bins, contigs) = gather_initial_bins_and_contig_information(&TEST_DATA_DIR, bin_gener);
        assert_eq!(bins.len(), 4);
        for bin in bins {
            match bin.bin_contigs {
                Some(the_contigs) => println!("{}", the_contigs.len()),
                None => panic!("Could not find contigs..")
            }
        }
        assert_eq!(contigs.len(), 1831);

    }
}