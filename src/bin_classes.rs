use std::{ffi::OsString, path::{PathBuf, Path}, fs::{File, self}, io::{Read, Write}, process::Command, sync::Arc, hash::{self, Hasher}, thread, time::Duration, collections::{HashSet, HashMap}};
use crate::{contigs::{Contig, EukaryoticContigInformation, ProkaryoticContigInformation}, utils::{get_contig_headers_from_fasta_file, generate_hash_from_contig_comb}, bin_info_storing::BinInfoStorer};
use itertools::Itertools;
use log::{debug, warn, info};
use std::hash::{Hash};
pub struct BinGenerator {
    pub eukaryotic_bin_creator: Option<CompleasmBinQualityGetter>,
    pub prokaryotic_bin_creator: Option<Checkm2BinQualityGetter>,
    pub hash_directory: PathBuf,
    pub maximum_contamination: f64,
    pub minimum_completeness: f64,


}

impl BinGenerator {
    pub fn initialise_bin_gen(checkm2_db_path: Option<PathBuf>, compleasm_lib_db_path: Option<PathBuf>, hash_directory: PathBuf, maximum_contamination: f64, minimum_completeness: f64) -> BinGenerator {

        let mut bin_gen = BinGenerator{
            eukaryotic_bin_creator: None,
            prokaryotic_bin_creator: None,
            hash_directory: hash_directory,
            maximum_contamination: maximum_contamination,
            minimum_completeness: minimum_completeness,

        };

    match checkm2_db_path {
        Some(path) => bin_gen.prokaryotic_bin_creator = Some(Checkm2BinQualityGetter { checkm2_db_path: path }),
        None => ()
        }
    match compleasm_lib_db_path {
        Some(path) => bin_gen.eukaryotic_bin_creator = Some(CompleasmBinQualityGetter { compleasm_db_lib_path: path, number_of_marker_ids: 255 }),
        None => ()
        }
    bin_gen
    }


    pub fn create_bin_from_original_bin_fasta_file(&self, fasta_path: &PathBuf) -> Option<Bin> {

        let contig_headers = get_contig_headers_from_fasta_file(fasta_path);
        let bin_hash = generate_hash_from_contig_comb(&contig_headers);
        let bin_directory_hash = &self.hash_directory.join(&bin_hash);
        println!("Original bin file path is: {}, hash is: {} ", fasta_path.to_str().unwrap(), &bin_hash);
        let new_fasta_path = &bin_directory_hash.join((&format!("{}.fa", bin_hash)));
        if !bin_directory_hash.is_dir() {
            fs::create_dir(&bin_directory_hash).unwrap();
            fs::copy(fasta_path, new_fasta_path).unwrap();
        } else {
            if !new_fasta_path.exists() {
                panic!("Found a premade bin directory with no fasta path in it. Likely the program crashed/was stopped and have faulty contig hash directories...");
            }
        }

        let bin_type = BinTypePredictor::run_eukrep_to_get_bin_type(&new_fasta_path, bin_directory_hash).unwrap();
        match self.get_bin_quality(&bin_type, &new_fasta_path, &bin_directory_hash, None) {
            Some(results) => Some(Bin {
                fasta_path: new_fasta_path.clone(),
                bin_contigs: None,
                completeness: Some(results.0),
                contamination: Some(results.1),
                bin_type: Some(bin_type),
                bin_hash: bin_hash.to_string(),
                bin_dir_path: bin_directory_hash.clone()
            }),
            None => None


        }


    }

    fn check_bin_reaches_minimum_quality(&self, comp_and_cont_values: &(f64, f64)) -> bool {
        if comp_and_cont_values.0 < self.minimum_completeness {
            return false
        }
        if comp_and_cont_values.1 > self.maximum_contamination {
            return false
        }
        true
    }

    pub fn create_bin_struct_from_contigs(&self, bin_contigs: Vec<Arc<Contig>>) -> Option<Bin> {
        let contig_headers = bin_contigs.iter().map(|contig| contig.header.clone()).collect_vec();
        let bin_hash = generate_hash_from_contig_comb(&contig_headers);
        let bin_directory_hash = &self.hash_directory.join(&bin_hash);
        let fasta_path = bin_directory_hash.join(&format!("{}.fa", &bin_hash));
        match self.get_bin_type_and_quality_for_constructed_bin(bin_directory_hash) {
            Some(qual_result) => {
                if self.check_bin_reaches_minimum_quality(&(qual_result.1, qual_result.2)) {
                
                Some(Bin {
                    fasta_path: fasta_path.clone(),
                    bin_contigs: Some(bin_contigs),
                    completeness: Some(qual_result.1),
                    contamination: Some(qual_result.2),
                    bin_type: Some(qual_result.0),
                    bin_hash: bin_hash,
                    bin_dir_path: bin_directory_hash.clone()
                })
            } else {
                None

                }
        
            },
            None => None
        }

    }

    fn get_bin_quality(&self, bin_type: &BinType, fasta_path: &PathBuf, bin_directory_hash: &PathBuf, contigs: Option<&Vec<Arc<Contig>>> ) -> Option<(f64, f64)> {
        let mut quality_result = None;
        match bin_type {
            BinType::eukaryote => {
                match &self.eukaryotic_bin_creator {
                    Some(euk_qual_getter) => {match euk_qual_getter.run_quality_check_for_bin(fasta_path, bin_directory_hash) {

                            Ok(result) => {
                            if self.check_bin_reaches_minimum_quality(&result) {

                                quality_result = Some((result.0, result.1))

                            } else {

                                quality_result = None
                            }

                            },

                            Err(_) => quality_result = None

                        }
                    }
                    None => quality_result = None

                }

            }

            BinType::prokaryote => {
                match &self.prokaryotic_bin_creator {
                    Some(prok_qual_getter) => {match prok_qual_getter.run_quality_check_for_bin(fasta_path, bin_directory_hash) {
                        Ok(result) => {
                            if self.check_bin_reaches_minimum_quality(&result) {
                                quality_result = Some((result.0, result.1))
                        } else {
                            quality_result = None
                        }

                        },
                        Err(_) => quality_result = None
                    }
                }
                    None => quality_result = None
                }

            }

        }
        let done_file_path = PathBuf::from(&bin_directory_hash.join("done_file.txt"));
        File::create(done_file_path).unwrap();
        quality_result

    }

    pub fn create_new_bin_from_contigs(&self, contigs: &Vec<Arc<Contig>>) -> Option<(String, f64, f64)> {
        let fasta_headers: Vec<String> = contigs.iter().map(|x| x.header.clone()).collect();
        let bin_hash = generate_hash_from_contig_comb(&fasta_headers);
        let bin_directory_hash = self.hash_directory.join(&bin_hash);
        let done_file_path = &bin_directory_hash.join("done_file.txt");
        let mut bin_info = None;
        if done_file_path.is_file() {
            println!("Done file path found!");
            return self.find_out_if_bin_has_already_been_produced_and_give_info(&bin_directory_hash, &bin_hash)

        }
        else {
            match fs::create_dir(&bin_directory_hash) {
                Ok(_) => return {

                    let new_bin_type = self.create_necessary_files_for_new_bin(contigs, &bin_directory_hash, &bin_hash);
                    let fasta_path = &bin_directory_hash.join(&format!("{}.fa", bin_hash));
                    let result = self.get_bin_quality(&new_bin_type, fasta_path, &bin_directory_hash, Some(contigs));
                    match result {

                        Some(qual_info) => Some((bin_hash, qual_info.0, qual_info.1)),
                        None => None

                    }

                },
                Err(_) => BinGenerator::wait_for_done_file(&bin_directory_hash)
            }
            bin_info = self.find_out_if_bin_has_already_been_produced_and_give_info(&bin_directory_hash, &bin_hash);
        }
        bin_info
    }

    fn find_out_if_bin_has_already_been_produced_and_give_info(&self, bin_dir_path: &PathBuf, bin_hash: &str) -> Option<(String, f64, f64)> {
        let done_file_path = bin_dir_path.join("done_file.txt");

        if done_file_path.is_file() {

            match self.get_bin_type_and_quality_for_constructed_bin(&bin_dir_path) {
                Some(bin_result) => {if self.check_bin_reaches_minimum_quality(&(bin_result.1, bin_result.2)) {
                    Some((bin_hash.to_string(), bin_result.1, bin_result.2))
                } else {
                        None
                    }
                },

                None => None
            }

        } else {
            None
        }
    }
    pub fn get_bin_type_and_quality_for_constructed_bin(&self, bin_dir_path: &PathBuf) -> Option<(BinType, f64, f64)> {
        let prok_qual_report_path = bin_dir_path.join("checkm2_results/quality_report.tsv");
        let euk_qual_report_path =  &bin_dir_path.join("compleasm_results/summary.txt");


        if prok_qual_report_path.is_file() {
            let quality = Checkm2BinQualityGetter::get_comp_and_cont_from_report_file_path(&prok_qual_report_path).unwrap();
            println!("Bin is prokaryote, completeness is: {}, contamination is: {}", quality.0, quality.1);
            Some((BinType::prokaryote, quality.0, quality.1))

        } else if euk_qual_report_path.is_file() {

            let quality = CompleasmBinQualityGetter::get_comp_and_cont_from_report_file_path(&euk_qual_report_path).unwrap();
            println!("Bin is eukaryote, completeness is: {}, contamination is: {}", quality.0, quality.1);
            Some((BinType::eukaryote, quality.0, quality.1))

        } else {
            None
        }
    }



    fn create_necessary_files_for_new_bin(&self, contigs: &Vec<Arc<Contig>>, hash_directory_path: &PathBuf, bin_hash: &str) -> BinType {

        let fasta_path = hash_directory_path.join(&format!("{}.fa", bin_hash));

        self.create_new_fasta_file(contigs, &fasta_path);

        let contig_prokaryotic_information: Vec<ProkaryoticContigInformation> = contigs.iter()
            .filter_map(|x| {match x.prokaryotic_contig_info.clone() {
                Some(prok_info) => Some(prok_info),
                None => None

            }}).collect();
        if (contigs.len() / 2) < contig_prokaryotic_information.len() { // ie if more than half of the contigs were designated as part of prokaryotic bins
            let diamond_lines_from_contigs: Vec<Vec<String>> = contig_prokaryotic_information.iter().map(|prok_info| prok_info.diamond_lines.clone()).collect();
            let protein_lines_from_contigs: Vec<Vec<String>> = contig_prokaryotic_information.iter().map(|prok_info| prok_info.predicted_protein_lines.clone()).collect();
            fs::create_dir(&hash_directory_path.join("checkm2_results/"));
            let diamond_dir = &hash_directory_path.join("checkm2_results/diamond_output/");
            fs::create_dir(diamond_dir);
            let diamond_path = &diamond_dir.join("DIAMOND_RESULTS.tsv");

            BinGenerator::create_new_diamond_file(diamond_lines_from_contigs, diamond_path, bin_hash);

            let protein_dir = &hash_directory_path.join("checkm2_results/protein_files/");
            fs::create_dir(protein_dir);
            let protein_path = &protein_dir.join(&format!("{}.faa", bin_hash));

            BinGenerator::create_new_protein_file(protein_lines_from_contigs, protein_path);
            return BinType::prokaryote
        }

        
        BinType::eukaryote


    }
    fn create_new_fasta_file(&self, contig_refs: &Vec<Arc<Contig>>, fasta_path: &PathBuf) {

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

    fn create_new_diamond_file(diamond_lines_from_contigs: Vec<Vec<String>>, diamond_path: &PathBuf, bin_hash: &str) {

        let diamond_vec: Vec<String> = diamond_lines_from_contigs.into_iter().flatten().map(|x| format!("{}Ω{}", &bin_hash, x)).collect();
        let diamond_string: String = diamond_vec.join("\n");
        let mut diamond_file = File::create(diamond_path).expect("can't create diamond file");
        diamond_file.write_all(diamond_string.as_bytes()).expect("Error could not write diamond string to file");

    }

    fn create_new_protein_file(protein_lines_from_contigs: Vec<Vec<String>>, protein_path: &PathBuf) {
        let protein_vec: Vec<String> = protein_lines_from_contigs.into_iter().flatten().collect();
        let protein_string: String = protein_vec.join("");
        let mut protein_file = File::create(protein_path).expect("can't create protein file");
        protein_file.write_all(protein_string.as_bytes()).expect("Error could not write protein string to file");
    }
    fn wait_for_done_file(bin_hash_dir: &PathBuf) {
        let done_file_path = &bin_hash_dir.join("done_file.txt");
        let mut time_out = 0;
        while !done_file_path.is_file() {
            time_out += 1;
            if time_out == 300000 {
                panic!("Contig hash directory for bin found, but waiting for done file still not generated. Panick. This was at: {}", &bin_hash_dir.to_string_lossy());
            }
            thread::sleep(Duration::from_millis(1));
        }
    }

}





#[derive(Clone, Debug)]
pub struct Bin {

    pub fasta_path: PathBuf,
    pub bin_contigs: Option<Vec<Arc<Contig>>>,
    pub completeness: Option<f64>,
    pub contamination: Option<f64>,
    pub bin_type: Option<BinType>,
    pub bin_hash: String,
    pub bin_dir_path: PathBuf


}

impl Bin {

 /*    fn add_contig_references_to_bin(&mut self, contigs: Vec<&'static Contig>) {
        self.bin_contigs = Some(contigs);
    } */

    pub fn add_contig_ref_to_bin(&mut self, contig: Arc<Contig>) {
        let bin_contigs = self.bin_contigs.clone();
        match bin_contigs {
            Some(current_contigs) => {
                let mut new_current_contigs = current_contigs.clone();
                new_current_contigs.push(contig);
                self.bin_contigs = Some(new_current_contigs.clone());
            }
            None => {
                self.bin_contigs = Some(vec![contig]);
            }
        }
    }
}

impl PartialEq for Bin {
    fn eq(&self, other: &Self) -> bool {
        self.bin_hash == other.bin_hash
    }
}

impl Eq for Bin {}

impl Hash for Bin {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.bin_hash.hash(state);
    }
}




struct BinTypePredictor;

impl BinTypePredictor {
    fn run_eukrep_to_get_bin_type(fasta_path: &PathBuf, output_directory: &PathBuf) -> Result<BinType, String> {
        let eukrep_res_path = PathBuf::from(output_directory.join("bin_eukrep_results.fa"));
        let eukrep_prok_res_path = PathBuf::from(output_directory.join("bin_prok_eukrep_results.fa"));

        if !eukrep_res_path.is_file() {


            let status = Command::new("EukRep")
            .args(["-i", fasta_path.to_str().unwrap(), "-o", eukrep_res_path.to_str().unwrap(), "--prokarya", eukrep_prok_res_path.to_str().unwrap()])
            .status()
            .map_err(|e| format!("Failed to run eukrep: {}", e))?;
            if status.success() {
                let bin_type_result = BinTypePredictor::give_bin_type_from_eukrep_result(fasta_path, &eukrep_res_path, &eukrep_prok_res_path).unwrap();
                Ok(bin_type_result)


            } else {
                Err(format!("Eukrep failed: {}", status))
                }
        } else {

            let bin_type_result = BinTypePredictor::give_bin_type_from_eukrep_result(fasta_path, &eukrep_res_path, &eukrep_prok_res_path).unwrap();
            Ok(bin_type_result)

        }
    }

    fn give_bin_type_from_eukrep_result(fasta_path: &PathBuf, eukrep_res_path: &PathBuf, eukrep_prok_res_path: &PathBuf) -> Result<BinType, String> {

        if let Ok(mut file) = File::open(eukrep_res_path) {

            if let Ok(mut prok_file) = File::open(eukrep_prok_res_path) {

                let mut eukaryotic_file_lines = String::new();
                file.read_to_string(&mut eukaryotic_file_lines).expect("Could not read eukrep file to string");
                let eukaryotic_header_count = eukaryotic_file_lines.lines().collect::<Vec<&str>>().len() / 2;
                let mut prok_file_lines = String::new();
                file.read_to_string(&mut prok_file_lines).expect("Could not read prok eukrep file to string");
                let prok_header_count = get_contig_headers_from_fasta_file(fasta_path).len() - eukaryotic_header_count;
                if prok_header_count > eukaryotic_header_count {
                    Ok(BinType::prokaryote)
                } else {
                    Ok(BinType::eukaryote)
                }

            } else {
                Err("could not open eukrep prokaryotic file".to_string())
            }
        } else {
            Err("Could not open eukaryotic file".to_string())
            }
    }


}




pub trait GetBinQuality {

    fn run_quality_check_for_bin(&self, fasta_path: &PathBuf, output_directory: &PathBuf) -> Result<(f64, f64), String>;
    fn get_comp_and_cont_from_report_file_path(quality_report_file_path: &PathBuf) -> Result<(f64, f64), String>;

}

pub struct CompleasmBinQualityGetter {
    compleasm_db_lib_path: PathBuf,
    number_of_marker_ids: usize
}

impl CompleasmBinQualityGetter {
    fn calculate_bin_quality_based_on_contig_buscos(&self, contigs: &Vec<Arc<Contig>>, output_file_path: &PathBuf) -> (f64, f64) {
    
        let contigs_with_ids: Vec<EukaryoticContigInformation> = contigs.iter().filter_map(|x| x.eukaryotic_contig_info.clone()).collect();
        let all_busco_ids: Vec<String> = contigs_with_ids.iter().map(|x| x.complete_buscos.clone()).flatten().collect();
        let number_of_ids = all_busco_ids.len();
        let unique_busco_ids: HashSet<String> = all_busco_ids.into_iter().collect();
        let perc_completeness = (unique_busco_ids.len() as f64 / self.number_of_marker_ids as f64) * 100.0;
        let number_of_duplicate_ids = number_of_ids - unique_busco_ids.len();
        let perc_contamination = (number_of_duplicate_ids as f64 / self.number_of_marker_ids as f64) * 100.0;

        self.write_brief_bin_qual_file(output_file_path, (perc_completeness, perc_contamination));
        (perc_completeness, perc_contamination)
    }
    
    fn write_brief_bin_qual_file(&self, output_file_path: &PathBuf, perc_comp_and_cont: (f64, f64)) {

        let mut output_file = File::create(output_file_path).unwrap();
        let string_for_file = format!("Completeness,{}\n Contamination,{}\ncalculated using {} markers", perc_comp_and_cont.0, perc_comp_and_cont.1, self.number_of_marker_ids);
        output_file.write_all(string_for_file.as_bytes()).unwrap(); 
    
    }

    fn get_qual_info_skipping_running_compleasm(&self, output_directory: &PathBuf, contigs: &Vec<Arc<Contig>>) -> Result<(f64, f64), String>{
        let quick_result_path = output_directory.join("compleasm_results/quick_summary.txt");
        if quick_result_path.is_file() {
            if let Ok(qual_result) = self.read_quick_summary_info(&quick_result_path) {
                Ok(qual_result)
            } else {
                Err(("Failed to get quick comp/cont".to_string()))
            }

        
        } else {
        
            Ok(self.calculate_bin_quality_based_on_contig_buscos(contigs, &quick_result_path))
        
        }

    }
    fn read_quick_summary_info(&self, quick_result_path: &PathBuf) -> Result<(f64, f64), String> {
        if let Ok(mut file) = File::open(quick_result_path){
            let mut qual_file_text = String::new();
            file.read_to_string(&mut qual_file_text).expect("can't read quick compleasm summary file");
            let file_lines: Vec<&str> = qual_file_text.lines().collect();
            let completeness = file_lines[0].split(",").last().unwrap();
            let contamination = file_lines[1].split(",").last().unwrap();
            Ok((completeness.parse::<f64>().unwrap(), contamination.parse::<f64>().unwrap()))

        }
        else {

            Err(("Could not find bin file when getting comp and contamination from quick summary compleasm report".to_string()))

            }
        }

    



}

impl GetBinQuality for CompleasmBinQualityGetter {
    fn run_quality_check_for_bin(&self, fasta_path: &PathBuf, output_directory: &PathBuf) -> Result<(f64, f64), String> {
        let bin_path_str = fasta_path.to_str().unwrap();
        let compleasm_directory = format!("{}/compleasm_results/", output_directory.to_string_lossy());
        let compleasm_result_path = output_directory.join("compleasm_results/summary.txt");
        if compleasm_result_path.is_file() {
            let (completeness, contamination) = Self::get_comp_and_cont_from_report_file_path(&compleasm_result_path).unwrap();
            return Ok((completeness, contamination))
        }
        let status = Command::new("compleasm")

            .args(["run", "-l", "eukaryota", "-L", &self.compleasm_db_lib_path.to_string_lossy(), "-a", bin_path_str, "-o", &compleasm_directory, "-t", "1", "-m", "lite"])
            .status()
            .map_err(|e| format!("Failed to execute compleasm: {}", e))?;

        if status.success() {
            let (completeness, contamination) = Self::get_comp_and_cont_from_report_file_path(&compleasm_result_path).unwrap();
            return Ok((completeness, contamination))
        }
        Err(format!("BAD: {}", status))
    }

    fn get_comp_and_cont_from_report_file_path(compleasm_file_path: &PathBuf) -> Result<(f64, f64), String> {
        if let Ok(mut file) = File::open(compleasm_file_path) {

            let mut compleasm_file_text = String::new();
            file.read_to_string(&mut compleasm_file_text).expect("Can't read compleasm file");
            let lines: Vec<&str> = compleasm_file_text.lines().collect();
            let completeness: f64 = lines[1].split("%").next().unwrap().split(":").skip(1).next().unwrap().parse().unwrap();
            let contamination: f64 = lines[2].split("%").next().unwrap().split(":").skip(1).next().unwrap().parse().unwrap();
            Ok((completeness, contamination))

        } else {

            Err("Failed to get completeness/contamination".to_string())

        }
    }

}

pub struct Checkm2BinQualityGetter {
    checkm2_db_path: PathBuf
}

impl GetBinQuality for Checkm2BinQualityGetter {

    fn run_quality_check_for_bin(&self, fasta_path: &PathBuf, output_directory: &PathBuf) -> Result<(f64, f64), String> {

        let bin_path_str = fasta_path.to_str().unwrap();
        let checkm2_analysis_result_directory = format!("{}/checkm2_results/", output_directory.to_string_lossy());
        let quality_report_path = PathBuf::from(&format!("{}/checkm2_results/quality_report.tsv", output_directory.to_string_lossy()));

        if quality_report_path.is_file() {
            let (completeness, contamination) = Self::get_comp_and_cont_from_report_file_path(&quality_report_path).unwrap();
            debug!("Found checkm2 quality report path!");
            return Ok((completeness, contamination))
        }
        let checkm2_db_path_str = &self.checkm2_db_path.to_string_lossy();
        let mut checkm2_args = vec!["predict", "-i", bin_path_str, "-x", "fa", "-o", &checkm2_analysis_result_directory, "--database_path", checkm2_db_path_str, "-t", "1", "--force", "--quiet"];
        if PathBuf::from(&checkm2_analysis_result_directory).is_dir() {
            checkm2_args.push("--resume");
        }
        else {
            fs::create_dir(&checkm2_analysis_result_directory);
        }
        let status = Command::new("checkm2")
            .args(checkm2_args)
            .status()
            .map_err(|e| format!("Failed to execute checkm2: {}", e))?;

        if status.success() {

            let (completeness, contamination) = Self::get_comp_and_cont_from_report_file_path(&quality_report_path).unwrap();
            Ok((completeness, contamination))

        } else {

            Err(format!("BAD: {}", status))

        }

    }

    fn get_comp_and_cont_from_report_file_path(quality_report_path: &PathBuf) -> Result<(f64, f64), String> {

        if let Ok(mut file) = File::open(quality_report_path){
            let mut qual_file_text = String::new();
            file.read_to_string(&mut qual_file_text).expect("can't read quality report file");
            let info_line = qual_file_text.lines().skip(1).next();
            match info_line {
                Some(qual_infos) => {
                    let completeness = qual_infos.split("\t").skip(1).next().unwrap().parse().unwrap();
                    let contamination = qual_infos.split("\t").skip(2).next().unwrap().parse().unwrap();
                    Ok((completeness, contamination))
                }
                None => Err(("Failure to find info line in checkm quality report!".to_string()))
            }

        }
        else {

            Err(("Could not find bin file when getting comp and contamination from quality report".to_string()))

            }
        }




}





#[cfg(test)]
mod tests {
    use crate::{initial_bins_and_contigs::BinContigInterface, contigs::ProkaryoticContigInformation};

    use super::*;
    use std::{fs, env, path::Path};
    use lazy_static::lazy_static;
    lazy_static! {
    static ref TEST_DATA_HASH: PathBuf = PathBuf::from("tests/test_data/bin_classes_testing/hash_results/");
}
    lazy_static! {
    static ref TEST_DATA_DIR: PathBuf = PathBuf::from("tests/test_data/bin_classes_testing/");
}

    lazy_static! { // unit test will only work if database is added here, not uploaded to git repo due to size
        static ref COMPLEASM_DB_LIB: PathBuf = PathBuf::from("tests/test_data/databases_for_testing/");
    }
    lazy_static! { // unit test will only work if database is added here, not uploaded to git repo due to size

        static ref CHECKM2_DB_PATH: PathBuf = PathBuf::from("tests/test_data/databases_for_testing/uniref100.KO.1.dmnd");

    }


    fn generate_bingen_struct() -> BinGenerator {
        let bin_gen = BinGenerator::initialise_bin_gen(
            Some(CHECKM2_DB_PATH.to_path_buf()), Some(COMPLEASM_DB_LIB.to_path_buf()), TEST_DATA_HASH.to_path_buf(), 0.0, 0.0);
            bin_gen
    }
    #[test]

    fn validate_prokaryote_bin_generation_and_premade() {
        let bin_gen = generate_bingen_struct();
        fs::remove_dir_all(&TEST_DATA_HASH.join("20ade716c5c524f277bb6a0cccc80e289494a5060df506f1090cf3736ac38b02/"));
        let bin_to_be_created = bin_gen.create_bin_from_original_bin_fasta_file(&PathBuf::from("tests/test_data/bin_classes_testing/premade_bin/Semibin_sample_1.fa")).unwrap();
        assert_eq!(bin_to_be_created.bin_type, Some(BinType::prokaryote));
        assert_eq!(bin_to_be_created.completeness, Some(100.0));
        assert_eq!(bin_to_be_created.contamination, Some(8.75));

        let premade_bin_to_be_created = bin_gen.create_bin_from_original_bin_fasta_file(&PathBuf::from("tests/test_data/bin_classes_testing/premade_bin/Semibin_sample_1.fa")).unwrap();
        assert_eq!(premade_bin_to_be_created.bin_type, Some(BinType::prokaryote));
        assert_eq!(premade_bin_to_be_created.completeness, Some(100.0));
        assert_eq!(premade_bin_to_be_created.contamination, Some(8.75));

    }

    fn validate_eukaryotic_bin_generation() {

        let bin_gen = generate_bingen_struct();
        let euk_test_fresh_bin = bin_gen.create_bin_from_original_bin_fasta_file(&PathBuf::from("tests/test_data/eukaryote_premade_bin/Saccharomyces_cerevisiae_eukaryote_test_sample.fa")).unwrap();
        assert_eq!(euk_test_fresh_bin.completeness, Some(96.74));
        assert_eq!(euk_test_fresh_bin.contamination, Some(0.00));

    }



    fn create_fake_prok_full_contigs_for_unit_test() -> Vec<Arc<Contig>> {
        // note these fakes are not consistent in their naming etc, more just to test the bin files are correctly being created
        let diamond_hypothetical_lines = vec!["ΩNODE_43646_length_1003_cov_2.925105_1\tUniRef100_W0U5L1~K00031\t81.7\t333\t60\t1\t1\t333\t23\t354\t4.4e-158\t562.4".to_string(),
        "ΩNODE_43757_length_1001_cov_2.635307_2\tUniRef100_A0A174C7U1~K04758\t78.6\t70\t15\t0\t1\t70\t1\t70\t1.5e-23\t113.2".to_string(),
        "ΩNODE_43443_length_1007_cov_2.412815_2\tUniRef100_D4IX81~K01649\t61.4\t83\t30\t2\t1\t83\t1\t81\t2.1e-20\t103.2".to_string()
        ];
        let protein_hypothetical_lines = vec![">NODE_43646_length_1003_cov_2.925105_1 # 660 # 872 # 1 # ID=1283_2;partial=00;start_type=ATG;rbs_motif=GGA/GAG/AGG;rbs_spacer=5-10bp;gc_cont=0.460\nMMPITMMKDGETGVIKKVGGKEETRKFLENLGFVVGGTVTVVSDIGGNLIVNVKDSRVAI\nPTKSARTIIAPDNMTGELQETALNVQYLKTVLG".to_string(),
        ">NODE_43646_length_1003_cov_2.925105_2 # 660 # 872 # 1 # ID=1283_2;partial=00;start_type=ATG;rbs_motif=GGA/GAG/AGG;rbs_spacer=5-10bp;gc_cont=0.460\nMMPITMMKDGETTCMKDDDEGVIKKVGGKEETRKFLENLGFVVGGTVTVVSDIGGNLIVNVKDSRVAI\nPTKSARTIIAPDNMTGELQETALNVQYLKTVLG".to_string(),
        ">NODE_43646_length_1003_cov_2.925105_3 # 660 # 872 # 1 # ID=1283_2;partial=00;start_type=ATG;rbs_motif=GGA/GAG/AGG;rbs_spacer=5-10bp;gc_cont=0.460\nMMPITMMKDGMMMDDDETTCMKDDDEGVIKKVGGKEETRKFLENLGFVVGGTVTVVSDIGGNLIVNVKDSRVAI\nPTKSARTIIAPDNMTGELQETALNVQYLKTVLG".to_string()
        ];
        let diamond_hypothetical_lines2 = vec!["ΩNODE_43647_length_1003_cov_2.925105_1\tUniRef100_W0U5L1~K00031\t81.7\t333\t60\t1\t1\t333\t23\t354\t4.4e-158\t562.4".to_string(),
        "ΩNODE_43647_length_1003_cov_2.925105_2\tUniRef100_A0A174C7U1~K04758\t78.6\t70\t15\t0\t1\t70\t1\t70\t1.5e-23\t113.2".to_string(),
        "ΩNODE_43647_length_1003_cov_2.925105_3\tUniRef100_D4IX81~K01649\t61.4\t83\t30\t2\t1\t83\t1\t81\t2.1e-20\t103.2".to_string()
        ];
        let protein_hypothetical_lines2 = vec![">NODE_43647_length_1003_cov_2.925105_1 # 660 # 872 # 1 # ID=1283_2;partial=00;start_type=ATG;rbs_motif=GGA/GAG/AGG;rbs_spacer=5-10bp;gc_cont=0.460\nMMPITDDDMMKDGETGVIKKVGGKEETRKFLENLGFVVGGTVTVVSDIGGNLIVNVKDSRVAI\nPTKSARTIIAPDNMTGELQETALNVQYLKTVLG".to_string(),
        ">NODE_43647_length_1003_cov_2.925105_2 # 660 # 872 # 1 # ID=1283_2;partial=00;start_type=ATG;rbs_motif=GGA/GAG/AGG;rbs_spacer=5-10bp;gc_cont=0.460\nMMPITMMKDGETTCMKDDDEGVIKDDKVGGKEETRKFLENLGFVVGGTVTVVSDIGGNLIVNVKDSRVAI\nPTKSARTIIAPDNMTGELQETALNVQYLKTVLG".to_string(),
        ">NODE_43647_length_1003_cov_2.925105_3 # 660 # 872 # 1 # ID=1283_2;partial=00;start_type=ATG;rbs_motif=GGA/GAG/AGG;rbs_spacer=5-10bp;gc_cont=0.460\nMMPITMMKDGMMMDDDETTCMDDDKDDDEGVIKKVGGKEETRKFLENLGFVVGGTVTVVSDIGGNLIVNVKDSRVAI\nPTKSARTIIAPDNMTGELQETALNVQYLKTVLG".to_string()
        ];
        let fake_contig_1 = Contig {
            header: ">NODE_43646_length_1003_cov_2.925105".to_string(),
            sequence: "ATTATATATAGTGCTAGCTACGATCGTAGCTACGTACGAGCTAGCATGCA".to_string(),
            prokaryotic_contig_info: Some(ProkaryoticContigInformation {diamond_lines: diamond_hypothetical_lines, predicted_protein_lines: protein_hypothetical_lines}),
            eukaryotic_contig_info: None,
            prok_or_euk: None
        };
        let fake_contig_2 = Contig {
            header: ">NODE_43647_length_1003_cov_2.925105".to_string(),
            sequence: "ATTATATATAGTGCTATTTATTATATAGCTACGATCGTAGCTACGTACGAGCTAGCATGCA".to_string(),
            prokaryotic_contig_info: Some(ProkaryoticContigInformation {diamond_lines: diamond_hypothetical_lines2, predicted_protein_lines: protein_hypothetical_lines2}),
            eukaryotic_contig_info: None,
            prok_or_euk: None
        };


        let contig_vec = vec![fake_contig_1, fake_contig_2];
        let fake_contig_arc: Vec<Arc<Contig>> = contig_vec.into_iter().map(|x| Arc::new(x)).collect();
        fake_contig_arc

    }

    /*
    #[test]
    fn validate_premade_bin() {

        let mut premade_bin_to_be_created = Bin {fasta_path: PathBuf::from("tests/test_data/hash_results/20ade716c5c524f277bb6a0cccc80e289494a5060df506f1090cf3736ac38b02.fa"),
            bin_contigs: None,
            completeness: None,
            contamination: None,
            bin_type: None,
            bin_hash: "20ade716c5c524f277bb6a0cccc80e289494a5060df506f1090cf3736ac38b02".to_string(),
            bin_dir_path: PathBuf::from("tests/test_data/hash_results/20ade716c5c524f277bb6a0cccc80e289494a5060df506f1090cf3736ac38b02/")

        };
        premade_bin_to_be_created.analyse_bin();
        assert_eq!(premade_bin_to_be_created.bin_type, Some(BinType::prokaryote));
        assert_eq!(premade_bin_to_be_created.completeness, Some(100.0));
        assert_eq!(premade_bin_to_be_created.contamination, Some(8.75));

    }
    #[test]
    fn validate_premade_bin_using_original_bin_file_to_bin_struct() {
        let bin_gener = create_bin_generator_for_testing();
        let test_premade_bin = bin_gener.create_bin_from_original_bin_fasta_file(&PathBuf::from("tests/test_data/bin_classes_testing/premade_bin/Semibin_sample_1.fa"));
             assert_eq!(test_premade_bin.bin_type, Some(BinType::prokaryote));
             assert_eq!(test_premade_bin.completeness, Some(100.0));
             assert_eq!(test_premade_bin.contamination, Some(8.75));

    }
    #[test]
    fn validate_fresh_bin_gen_process_using_original_bin_file_to_bin_struct() {
        let bin_gener = create_bin_generator_for_testing();

        let test_fresh_bin = bin_gener.create_bin_from_original_bin_fasta_file(&PathBuf::from("tests/test_data/premade_bin/Semibin_sample_1.fa"));
        assert_eq!(test_fresh_bin.bin_type, Some(BinType::prokaryote));
        assert_eq!(test_fresh_bin.completeness, Some(100.0));
        assert_eq!(test_fresh_bin.contamination, Some(8.75));


    }
    fn validate_eukaryote_premade_bin() {

        let bin_gener = create_bin_generator_for_testing();
        let euk_test_fresh_bin = bin_gener.create_bin_from_original_bin_fasta_file(&PathBuf::from("tests/test_data/eukaryote_premade_bin/Saccharomyces_cerevisiae_eukaryote_test_sample.fa"));
        assert_eq!(euk_test_fresh_bin.completeness, Some(96.74));
        assert_eq!(euk_test_fresh_bin.contamination, Some(0.00));

    }

    #[test]
    fn validate_create_new_bin_from_contigs() {
        let bin_gener = create_bin_generator_for_testing();
        let fake_contig_arcs = create_fake_prok_full_contigs_for_unit_test();
        let temp_test_dir = TEST_DATA_DIR.join("temp_bin_classes_test/");
        let fake_bin_hash_temp_path = &temp_test_dir.join("0a63bfe96821a441ef5e3f236716ccdd43d9072c95c60dbdb37a0982601555a1/");
        if temp_test_dir.is_dir() {
            fs::remove_dir(&temp_test_dir);
            fs::create_dir(&temp_test_dir);

        }
        fs::create_dir(fake_bin_hash_temp_path);
        bin_gener.create_necessary_files_for_new_bin(&fake_contig_arcs, fake_bin_hash_temp_path, "0a63bfe96821a441ef5e3f236716ccdd43d9072c95c60dbdb37a0982601555a1");
        let expected_fasta = fake_bin_hash_temp_path.join("0a63bfe96821a441ef5e3f236716ccdd43d9072c95c60dbdb37a0982601555a1.fa");
        let expected_protein = fake_bin_hash_temp_path.join("checkm2_results/protein_files/0a63bfe96821a441ef5e3f236716ccdd43d9072c95c60dbdb37a0982601555a1.faa");
        let expected_diamond = fake_bin_hash_temp_path.join("checkm2_results/diamond_output/DIAMOND_RESULTS.tsv");
        assert_eq!(true, expected_fasta.is_file());
        assert!(expected_fasta.is_file());
        assert!(expected_protein.is_file());
        assert!(expected_diamond.is_file());

       /*  let correct_protein_file_lines = read_file_to_string_for_validation(&PathBuf::from(&format!("{}fake_protein_file.faa", TEST_DATA_DIR.to_string_lossy())));
        let correct_diamond_file_lines = read_file_to_string_for_validation(&PathBuf::from(&format!("{}fake_diamond_results.tsv", TEST_DATA_DIR.to_string_lossy())));
        let correct_fasta_file_lines = read_file_to_string_for_validation(&PathBuf::from(&format!("{}fake_fasta.fa", TEST_DATA_DIR.to_string_lossy())));
        let expected_fasta_file_lines = read_file_to_string_for_validation(&expected_fasta);
        let expected_protein_file_lines = read_file_to_string_for_validation(&expected_protein);
        let expected_diamond_file_lines = read_file_to_string_for_validation(&expected_diamond);
       // assert_eq!(correct_diamond_file_lines, expected_diamond_file_lines);
        assert_eq!(correct_protein_file_lines.len(), expected_protein_file_lines.len());
        assert_eq!(correct_fasta_file_lines, expected_fasta_file_lines);
        */

    }
    fn create_premade_bin_1() -> Bin {
        // 1284 contigs
        let premade_bin = Bin {
            fasta_path: TEST_DATA_DIR.join("08c929dd57d754fe8e56f2befc8db4daed8bded44d54525e970e9be228345b9c/08c929dd57d754fe8e56f2befc8db4daed8bded44d54525e970e9be228345b9c.fa"),
            bin_contigs: None,
            completeness: Some(51.74),
            contamination: Some(13.2),
            bin_type: Some(BinType::prokaryote),
            bin_hash: "08c929dd57d754fe8e56f2befc8db4daed8bded44d54525e970e9be228345b9c".to_string(),
            bin_dir_path: TEST_DATA_DIR.join("08c929dd57d754fe8e56f2befc8db4daed8bded44d54525e970e9be228345b9c/")
        };
        premade_bin
    }
    #[test]
    fn validate_bin_creation() {
        let premade_bin_1 = create_premade_bin_1();
        let premade_bin_contigs: Vec<Arc<Contig>> = BinContigInterface::generate_contigs_from_bin(&premade_bin_1).into_iter().map(|x| Arc::new(x)).collect();
        let temp_test_dir = TEST_DATA_DIR.join("temp_bin_classes_test/");
        println!("{}", &temp_test_dir.to_string_lossy());
        let bin_gener = create_bin_generator_for_testing();
        let result = bin_gener.create_new_bin_from_contigs(&premade_bin_contigs);

    }
    fn read_file_to_string_for_validation(file_path: &PathBuf) -> String {
        let mut file = File::open(&file_path).expect(&format!("Could not open protein file for bin at path: {}", &file_path.to_string_lossy()));
        let mut lines = String::new();
        file.read_to_string(&mut lines).expect("Could not read protein file to string");
        lines

    }
    */
}
