use std::{path::PathBuf, fs::File, io::{Read, Write}, sync::Arc, process::Command, collections::{HashSet, HashMap}};
use glob::glob;
use itertools::Itertools;
use log::info;
use crate::{contigs::{Contig, EukaryoticContigInformation}};






pub struct EukaryoticBinQualityGetter {
    compleasm_db_lib_path: PathBuf,
    number_of_marker_ids: usize,
    compleasm_db: String,
}

impl EukaryoticBinQualityGetter {
    pub fn initialise(compleasm_db_lib_path: &str, number_of_markers: usize, compleasm_db: String) -> EukaryoticBinQualityGetter {
    
        EukaryoticBinQualityGetter { compleasm_db_lib_path: PathBuf::from(compleasm_db_lib_path), number_of_marker_ids: number_of_markers, compleasm_db: compleasm_db }

    }

    pub fn add_euk_info_to_contigs_using_compleasm(&self, contigs_path: &PathBuf, contigs: &mut [Contig], output_directory: &PathBuf, threads: usize) {

        let busco_original_table_path = output_directory.join(&format!("compleasm_results/{}/full_table_busco_format.tsv", &self.compleasm_db));
        if !busco_original_table_path.is_file() {
            info!("Initialising compleasm bin quality getter and running compleasm with {} db", &self.compleasm_db);
        
            self.run_compleasm_on_all_contigs(contigs_path, output_directory, threads).unwrap();

        } else {
            info!("Found busco table path for contigs, skipping running compleasm!");
        }

        let org_contig_busco_info = self.get_orf_contig_busco_info(busco_original_table_path);
        EukaryoticBinQualityGetter::add_busco_info_to_contigs(contigs, org_contig_busco_info);

        self.write_amount_of_buscos_in_contig_data(contigs, output_directory);
    
    }
    


    fn get_orf_contig_busco_info(&self, busco_original_table_path: PathBuf) -> Vec<(String, String)> {

        let mut busco_info_file = File::open(&busco_original_table_path).expect(&format!("Could not open busco table file at path: {}", busco_original_table_path.to_string_lossy()));
        let mut busco_lines = String::new();
        busco_info_file.read_to_string(&mut busco_lines).expect(&format!("Could not read busco table file to string"));
        let complete_or_duped_lines: Vec<&str> = busco_lines.lines().filter(|x| x.contains("Complete") || x.contains("Duplicate")).collect();
        let mut contig_names_and_ids = Vec::new();
        for line in complete_or_duped_lines {
            let mut split_lines = line.split("\t");
            let busco_id = split_lines.next().unwrap();
            split_lines.next();
            let contig_name = split_lines.next().unwrap();
            contig_names_and_ids.push((contig_name.to_string(), busco_id.to_string()));
        }
        contig_names_and_ids
         //   let mut busco_table_file = File::open(&protein_file_path).expect(&format!("Could not open protein file for bin at path: {}", &protein_file_path.to_string_lossy()));
    
    }

    fn add_busco_info_to_contigs(contigs: &mut [Contig], contig_names_and_ids: Vec<(String, String)>) {

        for contig in contigs {
            let mut busco_ids_for_contig = Vec::new();
            for (name, id) in &contig_names_and_ids {
                if contig.header.contains(name) {
                    busco_ids_for_contig.push(id.clone());
                }
            }
            if busco_ids_for_contig.len() != 0 {
                contig.add_eukaryotic_busco_info(busco_ids_for_contig);
            }
    
        }
    
    }

    fn write_amount_of_buscos_in_contig_data(&self, contigs: &[Contig], output_directory: &PathBuf) {
        let mut busco_marker_dict: HashMap<String, usize> = HashMap::new();
        for contig in contigs {
            match contig.eukaryotic_contig_info.clone() {
                Some(euk_info) => {
                    let busco_ids_in_contig = euk_info.complete_buscos;
                    let unique_busco_ids_in_contig: Vec<&String> = busco_ids_in_contig.iter().unique().collect();
                    if unique_busco_ids_in_contig.len() < busco_ids_in_contig.len() {
                        info!("Warning - contig has multiple duplicate unique buscos");
                    }
                    for busco_id in busco_ids_in_contig {

                    
                    if let Some(count) = busco_marker_dict.get_mut(&busco_id) {
                        *count = *count + 1;
                    } else {
                        busco_marker_dict.insert(busco_id.clone(), 1);
                    }
                }                 
                }
                None => ()
            }
        }
        let mut output_file = File::create(&output_directory.join("buscos_in_contigs.txt")).unwrap();
        let mut busco_in_contig_string = format!("Total buscos in chosen marker db: {}\n Total number of unique buscos in data: {}\nbusco_id,id_count\n", 
            &self.number_of_marker_ids, &busco_marker_dict.len());
        for (id, count) in busco_marker_dict {
            busco_in_contig_string.push_str(&format!("{},{}\n", &id, count));
        }
        output_file.write_all(busco_in_contig_string.as_bytes()).unwrap();


    }

    pub fn analyse_bin(&self, contigs: &[Arc<Contig>], output_hash_directory: &PathBuf) -> Option<(f64, f64)> {
    
        let quick_summary_file_path = output_hash_directory.join("eukaryotic_quick_summary.txt");
        if !quick_summary_file_path.is_file() {
            self.get_comp_and_cont_based_on_contig_markers(contigs, &quick_summary_file_path)
            
        } else {
            self.read_quick_summary_info(&quick_summary_file_path)
        }
    }

    fn get_comp_and_cont_based_on_contig_markers(&self, contigs: &[Arc<Contig>], quick_summary_file_path: &PathBuf) -> Option<(f64, f64)> {
        
        let contigs_with_ids: Vec<EukaryoticContigInformation> = contigs.iter().filter_map(|x| x.eukaryotic_contig_info.clone()).collect();
        let all_busco_ids: Vec<String> = contigs_with_ids.iter().map(|x| x.complete_buscos.clone()).flatten().collect();
        let number_of_ids = all_busco_ids.len();
        let unique_busco_ids: HashSet<String> = all_busco_ids.into_iter().collect();
        let perc_completeness = (unique_busco_ids.len() as f64 / self.number_of_marker_ids as f64)* 100.0;
        let number_of_duplicate_ids = number_of_ids - unique_busco_ids.len();
        let perc_contamination = (number_of_duplicate_ids as f64 / self.number_of_marker_ids as f64) * 100.0;
        let perc_cont_2_dps = (perc_contamination * 100.0).round() / 100.0;
        let perc_comp_2_dps = (perc_completeness * 100.0).round() / 100.0;
        
        self.write_brief_bin_qual_file(quick_summary_file_path, (perc_comp_2_dps, perc_cont_2_dps));
        
        if perc_completeness > 0.0 {
    
            Some((perc_comp_2_dps, perc_cont_2_dps))
      
        } else {
      
            None
     
        }
    }
    
    fn write_brief_bin_qual_file(&self, output_file_path: &PathBuf, perc_comp_and_cont: (f64, f64)) {

        let mut output_file = File::create(output_file_path).unwrap();
        let string_for_file = format!("Completeness:{}\n Contamination:{}\nCalculated using {} markers using {} database", 
            perc_comp_and_cont.0, perc_comp_and_cont.1, self.number_of_marker_ids, self.compleasm_db);
        output_file.write_all(string_for_file.as_bytes()).unwrap(); 
    
    }



    fn run_compleasm_on_all_contigs(&self, contigs_path: &PathBuf, output_directory: &PathBuf, threads: usize) -> Result<(), String> {
        
        let contig_path_str = contigs_path.to_str().unwrap();
        let compleasm_directory = format!("{}/compleasm_results/", output_directory.to_string_lossy());
        let taxon_from_db = &self.compleasm_db.split("_").next().unwrap(); // split taxon_odb to taxon eg eukaryota_odb10 -> eukaryota for running compleasm
        let output = Command::new("compleasm")
        .args([
            "run", "-l", taxon_from_db, "-L", &self.compleasm_db_lib_path.to_string_lossy(), "-t", &threads.to_string(), "-a", contig_path_str, "-o", &compleasm_directory])
        .output()
        .map_err(|e| format!("Failed to execute compleasm: {}", e))?;
    
        if output.status.success() {
            Ok(())
        } else {
            Err(format!(
                "Compleasm command exited with status: {}. Stderr: {}",
                output.status,
                String::from_utf8_lossy(&output.stderr)
            ))
        }
    
    }
    

    
    pub fn read_quick_summary_info(&self, quick_result_path: &PathBuf) -> Option<(f64, f64)> {
        
        let mut file = File::open(quick_result_path).unwrap();
        let mut qual_file_text = String::new();
        file.read_to_string(&mut qual_file_text).unwrap();
        let file_lines: Vec<&str> = qual_file_text.lines().collect();
        let completeness = file_lines[0].split(":").last().unwrap().parse::<f64>().unwrap();
        let contamination = file_lines[1].split(":").last().unwrap().parse::<f64>().unwrap();
        if completeness > 0.0 {
            
            Some((completeness, contamination))
   
        } else {
            None
        }

    }



}


struct CompleasmDB { // (busco datasets)

    number_of_markers: usize,
    compleasm_db_path: PathBuf,
    db_name: String

}

impl CompleasmDB {
    fn initialise(name_of_db: String, path_to_db: Option<PathBuf>) {

    }

    fn get_number_of_markers_from_database_folder(db_path: PathBuf) {


    }

}



#[cfg(test)]
mod tests {
    use crate::{contigs::{ProkaryoticContigInformation, ContigType}};

    use super::*;
    use std::{fs, env, path::Path};
    use lazy_static::lazy_static;
    lazy_static! {
    static ref TEST_DATA_DIR: PathBuf = PathBuf::from("tests/new_tests/test_bins/");
    static ref EUK_CONTIG_GATHER_OUTPUT_DIR: PathBuf = PathBuf::from("tests/new_tests/eukaryotic_contig_gatherer_test/");
}

    lazy_static! { // unit test will only work if database is added here, not uploaded to git repo due to size
        static ref COMPLEASM_DB_LIB: PathBuf = PathBuf::from("tests/new_tests/databases_for_testing/eukaryota_odb10");
    }

    fn create_example_compleasm_bq_getter_for_testing() -> EukaryoticBinQualityGetter {
        EukaryoticBinQualityGetter { compleasm_db_lib_path: COMPLEASM_DB_LIB.to_path_buf(), number_of_marker_ids: 255, compleasm_db: "eukaryota_odb10".to_string() }
    }


    fn create_example_eukaryotic_contigs() -> Vec<Arc<Contig>> {
        let mut contig_vec = Vec::new();
        contig_vec.push(Arc::new(Contig {
            header: "na".to_string(),
            sequence: "na".to_string(),
            prokaryotic_contig_info: None,
            eukaryotic_contig_info: Some(EukaryoticContigInformation {complete_buscos: vec!["779909at2759".to_string()]}),
            prok_or_euk: Some(ContigType::Eukaryote)
        }));
        contig_vec.push(Arc::new(Contig {
            header: "na".to_string(),
            sequence: "na".to_string(),
            prokaryotic_contig_info: None,
            eukaryotic_contig_info: Some(EukaryoticContigInformation {complete_buscos: vec!["290630at2759".to_string(), "290630at2759".to_string()]}), // duplicates so should only count once
            prok_or_euk: Some(ContigType::Eukaryote)
        }));


        contig_vec
    }

    fn create_example_eukaryotic_contigs_2() -> Vec<Arc<Contig>> {
        let mut contig_vec = Vec::new();
        contig_vec.push(Arc::new(Contig {
            header: "na".to_string(),
            sequence: "na".to_string(),
            prokaryotic_contig_info: None,
            eukaryotic_contig_info: Some(EukaryoticContigInformation {complete_buscos: vec!["779909at2759".to_string()]}),
            prok_or_euk: Some(ContigType::Eukaryote)
        }));
        contig_vec.push(Arc::new(Contig {
            header: "na".to_string(),
            sequence: "na".to_string(),
            prokaryotic_contig_info: None,
            eukaryotic_contig_info: Some(EukaryoticContigInformation {complete_buscos: vec!["290630at2759".to_string(), "290630at2759".to_string()]}), // duplicates so should only count once
            prok_or_euk: Some(ContigType::Eukaryote)
        }));


        contig_vec
    }



    fn initialise_test_euk_gatherer() -> EukaryoticBinQualityGetter {
        EukaryoticBinQualityGetter::initialise(&COMPLEASM_DB_LIB.to_str().unwrap(), 255, "eukaryota_odb10".to_string())
    }

    fn module_full_test() {
        let eukaryotic_bin_qual_getter = initialise_test_euk_gatherer();
        let expected_path = &EUK_CONTIG_GATHER_OUTPUT_DIR.join("compleasm_sacch_test_result_dir/compleasm_results/eukaryota_odb10/full_table_busco_format.tsv");
        if !expected_path.is_file() {
            test_run_compleasm_on_all_contigs()
        }


    }

    #[test]
    fn test_quick_summary_info() {
        let test_comp_bq_getter = create_example_compleasm_bq_getter_for_testing();
        let mut example_file_path = TEST_DATA_DIR.join("example_quick_summary_info_path.txt");
        let mut example_file = File::create(&example_file_path).unwrap();
        let example_file_string = format!("Completeness,90.0\n Contamination,5.0\ncalculated using 255 markers");
        example_file.write_all(example_file_string.as_bytes()).unwrap();
        let result = test_comp_bq_getter.read_quick_summary_info(&example_file_path).unwrap();
        assert_eq!(result.0, 90.0);
        assert_eq!(result.1, 5.0);


    }
    #[test]
    fn test_get_orf_contig_busco_info_method() {
        // 3 fragmented 14 missing 238 complete
        let eukaryotic_bin_qual_getter = initialise_test_euk_gatherer();
        let busco_original_table_path = &EUK_CONTIG_GATHER_OUTPUT_DIR.join("compleasm_sacch_test_result_dir/compleasm_results/eukaryota_odb10/full_table_busco_format.tsv");
        let test_complete_busco_ids = eukaryotic_bin_qual_getter.get_orf_contig_busco_info(busco_original_table_path.clone());
        if test_complete_busco_ids.len() != 238 {
            panic!("test_busco failed: Failed to get right number of busco ids")
        }
        let mut expected_id_file = File::open(&EUK_CONTIG_GATHER_OUTPUT_DIR.join("compleasm_sacch_test_result_dir/expected_ids.txt")).unwrap();
        let mut expected_id_string = String::new();
        expected_id_file.read_to_string(&mut expected_id_string);
        let expected_ids: Vec<String> = expected_id_string.lines().map(|line| line.to_string()).collect();

        for test_busco in test_complete_busco_ids {
            let busco_id = test_busco.1;
            if !expected_ids.contains(&busco_id) {
                panic!("test_busco failed: Test ids and expected ids don't match")
            }
        }

    }

    #[test]
    fn test_run_compleasm_on_all_contigs() {
        let compleasm_bin_getter = create_example_compleasm_bq_getter_for_testing();
        let contig_example_path = TEST_DATA_DIR.join("Saccharomyces_cerevisiae_eukaryote_test_sample.fa"); // example eukaryotic file
        let eukaryotic_bin_qual_getter = initialise_test_euk_gatherer();
        let result = eukaryotic_bin_qual_getter.run_compleasm_on_all_contigs(&contig_example_path, &EUK_CONTIG_GATHER_OUTPUT_DIR.join("compleasm_sacch_test_result_dir/"), 5).unwrap();
        let expected_path = &EUK_CONTIG_GATHER_OUTPUT_DIR.join("compleasm_sacch_test_result_dir/compleasm_results/eukaryota_odb10/full_table_busco_format.tsv");
        if !expected_path.is_file() {
            panic!("Compleasm failed to produce file in expected place");
        }
    }
    #[test]
    fn test_read_quick_summary_info() {
        let eukaryotic_bin_qual_getter = initialise_test_euk_gatherer();
        let test_file_path = &EUK_CONTIG_GATHER_OUTPUT_DIR.join("compleasm_sacch_test_result_dir/example_quick_summary_info.txt");
        match eukaryotic_bin_qual_getter.read_quick_summary_info(test_file_path) {
            Some(test_info) => {
                assert_eq!(75.0, test_info.0);
                assert_eq!(12.0, test_info.1);
            },
            None => panic!("test_read_quick_summary_info failed: method returned None")
        }

    }
    #[test]
    fn test_get_comp_and_cont_based_on_contig_markers_and_analyse_bin() {
        let test_contigs = create_example_eukaryotic_contigs();
        let eukaryotic_bin_qual_getter = initialise_test_euk_gatherer();
        let mut test_file_path = &EUK_CONTIG_GATHER_OUTPUT_DIR.join("compleasm_sacch_test_result_dir/temp_get_comp_and_cont_based_on_contig_markers.txt");
        fs::remove_file(&test_file_path);

        match eukaryotic_bin_qual_getter.get_comp_and_cont_based_on_contig_markers(&test_contigs, &test_file_path) {
            Some(quality) => {
                assert_eq!(0.78, quality.0);
                assert_eq!(0.39, quality.1);
            },
            None => panic!("Failed to get comp and cont based on contig markers!")
        }

        eukaryotic_bin_qual_getter.analyse_bin(&test_contigs, &EUK_CONTIG_GATHER_OUTPUT_DIR.join("compleasm_sacch_test_result_dir/"));
    }

}