use std::{path::PathBuf, process::Command, fs::File, io::Read};

use log::info;

use crate::{contigs::{Contig, ContigType}, utils::get_contig_headers_from_fasta_file};

pub struct ContigTypePredictor;


impl ContigTypePredictor {
    
    pub fn predict_contig_types(contigs_path: &PathBuf, contigs: &mut [Contig], output_directory: &PathBuf) -> Result<(usize, usize), String> {
        
        let eukrep_res_path = PathBuf::from(output_directory.join("all_contig_eukrep_results.fa"));
        let eukrep_prok_res_path = PathBuf::from(output_directory.join("all_contig_prok_eukrep_results.fa"));

        if !eukrep_res_path.is_file() || !eukrep_prok_res_path.is_file() {

 
            let status = Command::new("EukRep")
            .args(["-i", contigs_path.to_str().unwrap(), "-o", eukrep_res_path.to_str().unwrap(), "--prokarya", eukrep_prok_res_path.to_str().unwrap()])
            .status()
            .map_err(|e| format!("Failed to run eukrep: {}", e))?;
            if status.success() {
               
                match ContigTypePredictor::get_eukaryotic_and_prokaryotic_contig_headers(&eukrep_prok_res_path, &eukrep_res_path) {
                    
                    Ok((prokaryotic_headers, eukaryotic_headers)) => {
                    
                        ContigTypePredictor::add_type_to_contigs(contigs, &prokaryotic_headers, &eukaryotic_headers);
                        Ok((prokaryotic_headers.len(), eukaryotic_headers.len())) }
                    
                    Err(error_string) => Err(error_string)
                }

            } else {
                Err(format!("Eukrep failed: {}", status))
                }
        
        } else {
            
            match ContigTypePredictor::get_eukaryotic_and_prokaryotic_contig_headers(&eukrep_prok_res_path, &eukrep_res_path) {
                    
                Ok((prokaryotic_headers, eukaryotic_headers)) => {
                    
                    ContigTypePredictor::add_type_to_contigs(contigs, &prokaryotic_headers, &eukaryotic_headers);
                    Ok((prokaryotic_headers.len(), eukaryotic_headers.len())) }
                
                Err(error_string) => Err(error_string)
            
            }


        }
    }
    
    fn get_eukaryotic_and_prokaryotic_contig_headers(eukrep_prok_res_path: &PathBuf, eukrep_res_path: &PathBuf) -> Result<(Vec<String>, Vec<String>), String> {
        if eukrep_prok_res_path.is_file() && eukrep_res_path.is_file() {
        
            let prokaryotic_headers = get_contig_headers_from_fasta_file(eukrep_prok_res_path);
            let eukaryotic_headers = get_contig_headers_from_fasta_file(eukrep_res_path);
            info!("Eukrep identified {} contigs as prokaryote, and {} contigs as eukaryote", prokaryotic_headers.len(), eukaryotic_headers.len() );
            Ok((prokaryotic_headers, eukaryotic_headers)) 
        } else {
            Err("Could not find eukrep eukaryotic and prokaryotic files - likely something went wrong in eukrep!".to_string())
        }
    }



    fn add_type_to_contigs(contigs: &mut [Contig], prokaryotic_headers: &[String], eukaryotic_headers: &[String]) {
        for contig in contigs {
            if (eukaryotic_headers.contains(&contig.header) && prokaryotic_headers.contains(&contig.header)) || (!eukaryotic_headers.contains(&contig.header) && !prokaryotic_headers.contains(&contig.header)) { // the first one should never happen, but put just in case
                contig.prok_or_euk = None
            } else if eukaryotic_headers.contains(&contig.header) && !prokaryotic_headers.contains(&contig.header) {
                contig.prok_or_euk = Some(ContigType::Eukaryote)
            } else if !eukaryotic_headers.contains(&contig.header) && prokaryotic_headers.contains(&contig.header) {
                contig.prok_or_euk = Some(ContigType::Prokaryote)
            }

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
    static ref CONTIG_TYPE_PREDICTOR_DIR: PathBuf = PathBuf::from("tests/new_tests/contig_type_predictor_test/");
}


    #[test]
    fn test_contig_type_prediction() {
        
        
        let test_contigs_file_path = &TEST_DATA_DIR.join("all_unique_contigs.fa");
        let mut all_contigs = generate_all_contigs_from_contig_file(test_contigs_file_path);
        let test_eukrep_res_dir = &CONTIG_TYPE_PREDICTOR_DIR.join("test_all_contigs_results/");
        
      //  fs::remove_dir(&test_eukrep_res_dir).unwrap();
     //   fs::create_dir(test_eukrep_res_dir).unwrap();
        ContigTypePredictor::predict_contig_types(test_contigs_file_path, &mut all_contigs, test_eukrep_res_dir).unwrap();
        let contig_counts = all_contigs.iter().fold((0, 0, 0), |(euk, prok, unknown), contig| {
            match contig.prok_or_euk {
                Some(contig_type) => match contig_type {
                    ContigType::Eukaryote => (euk + 1, prok, unknown),
                    ContigType::Prokaryote => (euk, prok + 1, unknown),
                },
                None => (euk, prok, unknown + 1),
            }
        });
        println!("{:?}", contig_counts);
        assert_eq!(contig_counts.0, 30);
        assert_eq!(contig_counts.1, 1474)

    }



}

