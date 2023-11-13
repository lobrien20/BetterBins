use std::{path::PathBuf, process::Command, fs::{self, File}, io::{Read, Write}, sync::Arc, hash};

use log::{debug, info};

use crate::{contigs::{Contig, ProkaryoticContigInformation}, utils::create_new_fasta_file};



pub struct ProkaryoticBinQualityGetter {
    checkm2_db_path: PathBuf
}


impl ProkaryoticBinQualityGetter {
    pub fn initialise(checkm2_db_path_str: &str) -> ProkaryoticBinQualityGetter {
        ProkaryoticBinQualityGetter { checkm2_db_path: PathBuf::from(checkm2_db_path_str) }
    }

    pub fn add_prok_info_to_contigs_using_checkm2(&self, contigs_path: &PathBuf, contig_res_output_directory: &PathBuf, all_contigs: &mut [Contig], threads: usize) {
        
        let diamond_and_protein_paths = ProkaryoticBinQualityGetter::run_quality_check_for_all_contigs_to_gather_info(contigs_path, contig_res_output_directory, 
            self.checkm2_db_path.to_str().unwrap(), threads).unwrap();
        
        ProkaryoticBinQualityGetter::add_prokaryotic_orf_info_to_contigs(all_contigs, diamond_and_protein_paths);
    }
    
    fn run_quality_check_for_all_contigs_to_gather_info(contigs_path: &PathBuf, output_directory: &PathBuf, checkm2_db_path_str: &str, threads: usize) -> Result<(PathBuf, PathBuf), String> {

        // Technically this part can be skipped to just run prodigal/diamond since no need to run the checkm2 contigs. This is just a more simple way whilst also being safer in case
        // checkm2 gets updated/changed at all to have different settings for the model preparation steps

        let contig_path_str = contigs_path.to_str().unwrap();
        let checkm2_contig_result_directory = format!("{}/checkm2_results/", output_directory.to_string_lossy());
        match ProkaryoticBinQualityGetter::find_diamond_and_protein_files(&checkm2_contig_result_directory, contig_path_str) {
            Some(protein_and_diamond_paths) => {info!("Found pre-existing contig diamond and protein info"); return Ok(protein_and_diamond_paths)},
            None => (),
        }
        let thread_str = threads.to_string();
        fs::remove_dir(&checkm2_contig_result_directory);
        fs::create_dir(&checkm2_contig_result_directory).unwrap();
        let mut checkm2_args = vec!["predict", "-i", contig_path_str, "-x", "fa", "-o", &checkm2_contig_result_directory, "--database_path", checkm2_db_path_str, "-t", &thread_str, "--force", "--quiet"];

        info!("Running checkm with {} threads", &thread_str);
        let status = Command::new("checkm2")
            .args(checkm2_args)
            .status()
            .map_err(|e| format!("Failed to execute checkm2: {}", e))?;

        if status.success() {

            match ProkaryoticBinQualityGetter::find_diamond_and_protein_files(&checkm2_contig_result_directory, contig_path_str) {
                Some(protein_and_diamond_paths) => Ok(protein_and_diamond_paths),
                None => Err("Checkm2 of all contigs complete but couldn't find diamond and protein files!".to_string())
            }

        } else {

            Err(format!("BAD: {}", status))

        }


    }

    fn find_diamond_and_protein_files(checkm2_contig_result_directory: &str, contig_path_str: &str) -> Option<(PathBuf, PathBuf)> {
        
        let expected_diamond_info_path = PathBuf::from(&checkm2_contig_result_directory).join("diamond_output/DIAMOND_RESULTS.tsv");
        let contig_file_name = contig_path_str.split("/").last().unwrap();
        let split: Vec<&str> = contig_file_name.rsplitn(2, '.').collect();
        let contig_file_name_without_ext = split.last().unwrap_or(&"");
        println!("contig file name is: {}", contig_file_name_without_ext);
        let expected_protein_info_path = PathBuf::from(&checkm2_contig_result_directory).join(&format!("protein_files/{}.faa", contig_file_name_without_ext));
        println!("Expected diamond path path is: {}", &expected_diamond_info_path.to_str().unwrap());
        if expected_diamond_info_path.is_file() && expected_protein_info_path.is_file() {
            Some((expected_diamond_info_path, expected_protein_info_path))
        } else {
            None
        }
    }

    fn add_prokaryotic_orf_info_to_contigs(all_contigs: &mut [Contig], diamond_and_protein_paths: (PathBuf, PathBuf) ) {
        
        let protein_info = ProkaryoticBinQualityGetter::get_protein_info_from_file(&diamond_and_protein_paths.1);
        let diamond_info = ProkaryoticBinQualityGetter::get_prokaryotic_diamond_info_from_file(&diamond_and_protein_paths.0);
        let mut total_used_contigs_len = 0;
        for contig in all_contigs {

            let mut header_without_shark_symbol = contig.header.clone();
            header_without_shark_symbol.remove(0);
            let contig_protein_info: Vec<String> = protein_info.iter().filter(|x| x.contains(&contig.header)).map(|x| x.to_string()).collect();
            total_used_contigs_len += contig_protein_info.len();
            let contig_diamond_info: Vec<String> = diamond_info.iter().filter(|x| x.contains(&header_without_shark_symbol)).map(|x| x.to_string()).collect();
            if contig_protein_info.len() != 0 || contig_diamond_info.len() != 0 {
                contig.add_diamond_and_protein_information(contig_diamond_info, contig_protein_info);
            }
        }
        println!("total used protein cont: {}", total_used_contigs_len);

    }


    fn get_protein_info_from_file( protein_file_path: &PathBuf) -> Vec<String> {
    
        let mut protein_file = File::open(&protein_file_path).expect(&format!("Could not open protein file for bin at path: {}", &protein_file_path.to_string_lossy()));
        let mut protein_lines = String::new();
        protein_file.read_to_string(&mut protein_lines).expect("Could not read protein file to string");
        let mut protein_vec: Vec<String> = protein_lines.split(">").map(|x| x.to_string()).collect();
        println!("1 is: {}", protein_vec[0]);
        println!("2 is: {}", protein_vec[1]);
        println!("There are {}", protein_vec.len());
        

        let protein_info: Vec<String> = protein_vec.iter().skip(1).map(|x| format!(">{}", x)).collect();
        println!("There are {}", protein_info.len());
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

    pub fn analyse_bin(&self, contigs: &[Arc<Contig>], fasta_path: &PathBuf, hash_directory_path: &PathBuf, bin_hash: &str) -> Option<(f64, f64)> {
        println!("Analysing prokaryotic bin!");
        let checkm2_analysis_result_directory = hash_directory_path.join("checkm2_results/");
        if checkm2_analysis_result_directory.is_dir() {
        
            return self.get_comp_and_cont_from_report_file_path(&checkm2_analysis_result_directory.join("quality_report.tsv"))
        
        }

        fs::create_dir(&checkm2_analysis_result_directory).unwrap();
        
        match self.create_prokaryotic_specific_bin_files(contigs, &checkm2_analysis_result_directory, bin_hash) {
            Ok(_) => (),
            Err(_) => return None
        }
        let bin_path_str = fasta_path.to_str().unwrap();
        
        match self.run_checkm2_on_bin(bin_path_str, &checkm2_analysis_result_directory) {
            Some(quality) => Some(quality),
            None => None
        }

    }
    fn run_checkm2_on_bin(&self, bin_path_str: &str, checkm2_analysis_result_directory: &PathBuf) -> Option<(f64, f64)>  {


        let checkm2_db_path_str = &self.checkm2_db_path.to_string_lossy();
        let checkm2_args = vec!["predict", "-i", bin_path_str, "-x", "fa", "-o", &checkm2_analysis_result_directory.to_str().unwrap(), 
            "--database_path", checkm2_db_path_str, "-t", "1", "--force", "--quiet", "--resume"];

        let status = Command::new("checkm2")
            .args(checkm2_args)
            .status()
            .expect("Failed to execute checkm2");

        if status.success() {
            self.get_comp_and_cont_from_report_file_path(&PathBuf::from(&checkm2_analysis_result_directory).join("quality_report.tsv"))

        } else {
            
            info!("Checkm2 failure: {}", status);
            None

        }

    }


    fn create_prokaryotic_specific_bin_files(&self, contigs: &[Arc<Contig>], checkm2_analysis_result_directory: &PathBuf, bin_hash: &str) -> Result<(), String>  {

        let contig_prokaryotic_information: Vec<ProkaryoticContigInformation> = contigs.iter()
            .filter_map(|x| x.prokaryotic_contig_info.clone())
            .collect();
        
        if contig_prokaryotic_information.len() == 0 {
            return Err("No prokaryotic contig information".to_string())
        }
        let diamond_lines_from_contigs: Vec<Vec<String>> = contig_prokaryotic_information.iter().map(|prok_info| prok_info.diamond_lines.clone()).collect();
        let protein_lines_from_contigs: Vec<Vec<String>> = contig_prokaryotic_information.iter().map(|prok_info| prok_info.predicted_protein_lines.clone()).collect();
        if diamond_lines_from_contigs.len() == 1 && diamond_lines_from_contigs[0].len() == 0 {
            return Err("No prokaryotic diamonds".to_string())
        }
        println!("{:?}", diamond_lines_from_contigs);


        self.create_new_diamond_file(diamond_lines_from_contigs, &checkm2_analysis_result_directory, bin_hash);
        self.create_new_protein_file(protein_lines_from_contigs, &checkm2_analysis_result_directory, bin_hash);
        Ok(())

        


    }


    fn create_new_diamond_file(&self, diamond_lines_from_contigs: Vec<Vec<String>>, checkm2_analysis_result_directory: &PathBuf, bin_hash: &str) {

        let diamond_dir = &checkm2_analysis_result_directory.join("diamond_output/");
        fs::create_dir(diamond_dir).unwrap();
        let diamond_path = &diamond_dir.join("DIAMOND_RESULTS.tsv");

        let diamond_vec: Vec<String> = diamond_lines_from_contigs.into_iter().flatten().map(|x| format!("{}Ω{}", &bin_hash, x)).collect();
        let diamond_string: String = diamond_vec.join("\n");
        let mut diamond_file = File::create(diamond_path).expect("can't create diamond file");
        diamond_file.write_all(diamond_string.as_bytes()).expect("Error could not write diamond string to file");

    }

    fn create_new_protein_file(&self, protein_lines_from_contigs: Vec<Vec<String>>, checkm2_analysis_result_directory: &PathBuf, bin_hash: &str) {
        
        let protein_dir = &checkm2_analysis_result_directory.join("protein_files/");
        fs::create_dir(protein_dir).unwrap();
        let protein_path = &protein_dir.join(&format!("{}.faa", bin_hash));

        let protein_vec: Vec<String> = protein_lines_from_contigs.into_iter().flatten().collect();
        let protein_string: String = protein_vec.join("");
        let mut protein_file = File::create(protein_path).expect("can't create protein file");
        protein_file.write_all(protein_string.as_bytes()).expect("Error could not write protein string to file");
    
    }



    fn get_comp_and_cont_from_report_file_path(&self, quality_report_path: &PathBuf) -> Option<(f64, f64)> {

        if let Ok(mut file) = File::open(quality_report_path){
            let mut qual_file_text = String::new();
            file.read_to_string(&mut qual_file_text).expect("can't read quality report file");
            let info_line = qual_file_text.lines().skip(1).next();
            match info_line {
                Some(qual_infos) => {
                    let completeness = qual_infos.split("\t").skip(1).next().unwrap().parse().unwrap();
                    let contamination = qual_infos.split("\t").skip(2).next().unwrap().parse().unwrap();
                    Some((completeness, contamination))
                }
                None => None
            }

        }
        else {

            info!("Could not find bin file when getting comp and contamination from quality report");
            None

            }
        }


}


#[cfg(test)]
mod tests {
    use crate::{contigs::{ProkaryoticContigInformation, ContigType}, initialise_bins_and_contigs::generate_all_contigs_from_contig_file};

    use super::*;
    use std::{fs, env, path::Path};
    use itertools::Itertools;
    use lazy_static::lazy_static;
    use regex::Regex;
    lazy_static! {
    static ref TEST_DATA_DIR: PathBuf = PathBuf::from("tests/new_tests/test_bins/");
    static ref PROK_CONTIG_GATHERER_DIR: PathBuf = PathBuf::from("tests/new_tests/prokaryotic_contig_gatherer_test/");
}

    lazy_static! { // unit test will only work if database is added here, not uploaded to git repo due to size
        static ref CHECKM2_DB_PATH: PathBuf = PathBuf::from("tests/new_tests/databases_for_testing/uniref100.KO.1.dmnd");
    }


    #[test]
    fn test_add_prok_info_to_contigs() { // uses an example bin for this
        let prokaryotic_bin_quality_getter = ProkaryoticBinQualityGetter { checkm2_db_path: CHECKM2_DB_PATH.to_path_buf() };
        let test_contig_file_path = &TEST_DATA_DIR.join("Cow_36641_day_56_trimmed.066.fa");
        let mut test_bin_contigs = generate_all_contigs_from_contig_file(test_contig_file_path);
        println!("number of bin contigs: {}", test_bin_contigs.len());
        let checkm2_test_path = &PROK_CONTIG_GATHERER_DIR.join("checkm2_contigs_test/");
        prokaryotic_bin_quality_getter.add_prok_info_to_contigs_using_checkm2(test_contig_file_path, checkm2_test_path, &mut test_bin_contigs, 6);
        let mut total_proteins_predicted = 0;
        let mut total_diamonds_predicted = 0;
        for contig in test_bin_contigs {
            match contig.prokaryotic_contig_info {
                Some(prok_info) => {
                    let proteins_predicted = prok_info.predicted_protein_lines.iter()
                    .filter(|line| line.contains(">")).collect_vec();
                    total_proteins_predicted += proteins_predicted.len();
                    println!("{}", total_proteins_predicted);
                    total_diamonds_predicted += get_unique_diamonds_for_test(&prok_info).len();

                },
                None => ()
            }
        }
        assert_eq!(total_proteins_predicted, 2739);
        assert_eq!(total_diamonds_predicted, 689);
        //689 diamonds
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



    #[test]
    fn verify_get_comp_and_cont_from_report_file_path() {
        let prokaryotic_bin_quality_getter = ProkaryoticBinQualityGetter { checkm2_db_path: CHECKM2_DB_PATH.to_path_buf() };
        let test_file_path = &PROK_CONTIG_GATHERER_DIR.join("test_quality_report.tsv");
        match prokaryotic_bin_quality_getter.get_comp_and_cont_from_report_file_path(test_file_path) {
            Some(quality) => {
                assert_eq!(51.74, quality.0);
                assert_eq!(13.2, quality.1);
            },
            None => panic!("verify_getting_comp_and_cont_from_report_file_path test fail: Could not extract quality from report path")
        }
    }

}