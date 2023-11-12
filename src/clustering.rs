use std::{sync::Arc, collections::HashMap, path::PathBuf, fs::{File, self}, io::{Write, Read}, process::Command};

use blake3::Hash;
use itertools::Itertools;
use log::info;
use rayon::prelude::{IntoParallelRefIterator, ParallelIterator};

use crate::{bin_classes::{Bin, BinGenerator, BinType}, contigs::Contig};

// fs::create_dir(cluster_directory).unwrap_or(println!("Cluster results directory already exists."));

pub fn run_bin_clustering_module_to_generate_all_potential_bins(results_directory: &PathBuf, bins: Vec<Bin>, contigs: Vec<Arc<Contig>>, bin_generator: Arc<BinGenerator>) -> Vec<(Vec<Arc<Contig>>, String, f64, f64)> {
    let cluster_directory = results_directory.join("cluster_output/");
    fs::create_dir(&cluster_directory).unwrap_or(println!("Cluster results directory already exists."));
    let contig_presence_info_path = cluster_directory.join("contig_presence_info.tsv");
    let unique_bins = PrepareForClustering::remove_duplicate_bins(bins);
    if !contig_presence_info_path.is_file() {
      
        let contig_presence_hashmap = PrepareForClustering::create_contig_presence_hashmap(&unique_bins, &contigs);
        PrepareForClustering::write_contig_presence_info_to_file(&contig_presence_info_path, contigs.len(), contig_presence_hashmap)
    
    }
    let bin_clustered_hashmap = Clustering::cluster_bins(&contig_presence_info_path, &cluster_directory);
    info!("Creating and testing new hypothetical bins...");
    
    let bin_results = create_and_test_new_hypothetical_bins(bin_clustered_hashmap, &unique_bins, bin_generator, &cluster_directory);
    info!("{} New hypothetical bins created and tested - clustering stage complete!", bin_results.len());
    bin_results
    
}
fn create_and_test_new_hypothetical_bins(bin_clustered_hashmap: HashMap<String, usize>, bins: &Vec<Bin>, bin_generator: Arc<BinGenerator>, cluster_directory: &PathBuf) -> Vec<(Vec<Arc<Contig>>, String, f64, f64)> {
    
    let cluster_id_to_bin_hashmap = Clustering::generate_cluster_id_to_clustered_bins_hashmap(bin_clustered_hashmap, bins);
    
    let eukaryote_and_prokaryote_split_cluster_id_to_bin_hashmap = Clustering::filter_eukaryote_and_prokaryote_clusters(cluster_id_to_bin_hashmap);
    let filtered_clusters_output_file = cluster_directory.join("output_cluster_file_eukaryotes_and_prokaryotes_split.tsv");
    Clustering::create_new_cluster_file_from_eukaryote_and_prokaryote_filtering(&filtered_clusters_output_file, &eukaryote_and_prokaryote_split_cluster_id_to_bin_hashmap);

    Clustering::generate_viable_contig_combinations_for_each_cluster(eukaryote_and_prokaryote_split_cluster_id_to_bin_hashmap, bin_generator)
    
}
struct PrepareForClustering;


impl PrepareForClustering {
    fn remove_duplicate_bins(bins: Vec<Bin>) -> Vec<Bin> {
        let mut unique_bin_hashes = Vec::new();
        let mut unique_bins = Vec::new();
        for bin in bins {
            if !unique_bin_hashes.contains(&bin.bin_hash) {
                unique_bin_hashes.push(bin.bin_hash.clone());
                unique_bins.push(bin);
            
            } else {
                info!("Clustering: Found duplicate initial bin!");
            }
        }
        unique_bins
    }
    fn create_contig_presence_hashmap(bins: &Vec<Bin>, contigs: &Vec<Arc<Contig>>) -> HashMap<String, Vec<usize>>{
        let mut contig_presence_hashmap = HashMap::new();
        
        for bin in bins {
            match &bin.bin_contigs {
                Some(the_bin_contigs) => {
                    let contig_presence_vector = PrepareForClustering::create_contig_presence_vector_for_bin(the_bin_contigs, &contigs);
                    contig_presence_hashmap.insert(bin.bin_hash.clone(), contig_presence_vector);
                }
                None => {
                    panic!("Bin has no contigs. Error.")
                }
            }
        }
        println!("{:?}", contig_presence_hashmap);
        contig_presence_hashmap
    }

    fn create_contig_presence_vector_for_bin(bin_contigs: &Vec<Arc<Contig>>, contigs: &Vec<Arc<Contig>>) -> Vec<usize> {
        let mut contig_presence_vector = Vec::new();
        for contig in contigs {
            if bin_contigs.contains(contig) {
                contig_presence_vector.push(1);
            } else {
                contig_presence_vector.push(0);
            }
        }
        contig_presence_vector
    
    }
    fn write_contig_presence_info_to_file(output_path: &PathBuf, num_of_contigs: usize, contig_presence_hashmap: HashMap<String, Vec<usize>>) {
        let mut file_string = "bin_number".to_string();
        for i in 0..num_of_contigs {
            file_string.push_str(&format!("\tcontig_{i}").to_string());
        }
        file_string.push_str("\n");
        for (bin_num, presence_vals) in contig_presence_hashmap {
            let mut line_string = (&format!("{}", bin_num)).to_string();
            for val in presence_vals {
                line_string.push_str(&format!("\t{}", val));
            }
            line_string.push_str("\n");
            file_string.push_str(&line_string);
        }

        let mut contig_presence_info_file = File::create(output_path).expect("Can't create contig presence info file");

        contig_presence_info_file.write_all(file_string.as_bytes()).expect("Error could not write summariser string into summariser file");
    }
}

struct Clustering;

impl Clustering {
    fn cluster_bins(bin_contig_presence_file_path: &PathBuf, output_cluster_dir_path: &PathBuf) -> HashMap<String, usize> {
        let output_cluster_file_path = output_cluster_dir_path.join("output_cluster_file_path.tsv");
        if !bin_contig_presence_file_path.is_file() {
            panic!("Error, Cannot find cluster file path."); // fix this to become an error instead
        }
        if !output_cluster_file_path.is_file() {
            match Clustering::run_python_clustering_script(bin_contig_presence_file_path, &output_cluster_dir_path) {
                Some(cluster_hashmap) => return cluster_hashmap,
                None => panic!("Error, no cluster hashmap generated")
            }
        } else {
            Clustering::get_cluster_file_result(&output_cluster_file_path)
        }

    }

    fn run_python_clustering_script(input_contig_presences_path: &PathBuf, output_clusters_dir_path: &PathBuf) -> Option<HashMap<String, usize>> {
    
        let output = Command::new("python3")
        .args(["clustering.py", &input_contig_presences_path.to_string_lossy(), &output_clusters_dir_path.to_string_lossy()])
        .status()
        .map_err(|e| format!("Failed to run python clustering script: {}", e)).expect("Failed command");
        if output.success() {
            let bin_info_hashmap = Clustering::get_cluster_file_result(&output_clusters_dir_path.join("output_cluster_file_path.tsv"));
            Some(bin_info_hashmap)
        } else {
            None
        }
    }
    fn get_cluster_file_result(file_path: &PathBuf) -> HashMap<String, usize> {


        let mut file = File::open(file_path).unwrap();
        let mut file_text = String::new();
        file.read_to_string(&mut file_text).expect("Could not read file");
        let mut bin_info_hashmap: HashMap<String,usize>  = HashMap::new();
    
        for line in file_text.lines() {
            if line.contains("bin_id\tcluster_id") {
                continue
            }
            let bin_cluster_info: Vec<&str> = line.split("\t").collect();
            println!("{:?}", bin_cluster_info);
            bin_info_hashmap.insert(bin_cluster_info[0].parse::<String>().unwrap(), bin_cluster_info[1].parse::<usize>().unwrap());
        }
        bin_info_hashmap
    }
    fn generate_cluster_id_to_clustered_bins_hashmap<'a> (bin_info_hashmap: HashMap<String, usize>, bins: &'a Vec<Bin>) -> HashMap<usize, Vec<&'a Bin>>{
        let mut cluster_to_bin_hashmap: HashMap<usize, Vec<&Bin>> = HashMap::new();
        
        for (bin_hash, cluster_id) in bin_info_hashmap {

            for bin in bins {
                
                if bin.bin_hash == bin_hash {
                    
                    if let Some(bins) = cluster_to_bin_hashmap.get_mut(&cluster_id) {
                        bins.push(&bin);

                    } else {
                        let mut bins_vec = vec![bin];
                        cluster_to_bin_hashmap.insert(cluster_id, bins_vec);
                    }
                }
            }
        }
        cluster_to_bin_hashmap
    }

    fn generate_viable_contig_combinations_for_each_cluster(cluster_bin_hashmap: HashMap<usize, Vec<&Bin>>, bin_generator: Arc<BinGenerator>) -> Vec<(Vec<Arc<Contig>>, String, f64, f64)> {
        


        let all_contig_combination_results: Vec<Vec<(Vec<Arc<Contig>>, String, f64, f64)> > = cluster_bin_hashmap.values().cloned().collect::<Vec<Vec<&Bin>>>().par_iter().map(|bins_in_cluster| {
            let bin_gen_arc_ref = Arc::clone(&bin_generator);
            Clustering::generate_all_viable_bin_contig_combinations_for_cluster(bins_in_cluster, bin_gen_arc_ref)

        }).collect();
        
        let all_viable_contig_combinations: Vec<(Vec<Arc<Contig>>, String, f64, f64)>  = all_contig_combination_results.into_iter().flatten().collect();
        all_viable_contig_combinations

    }
    
    fn generate_all_viable_bin_contig_combinations_for_cluster(bins_in_cluster: &Vec<&Bin>, bin_gen_arc_ref: Arc<BinGenerator>) -> Vec<(Vec<Arc<Contig>>, String, f64, f64)> {
        let mut all_combinations: Vec<Vec<&Bin>> = (1..=bins_in_cluster.len()).flat_map(|i| bins_in_cluster.clone().into_iter().combinations(i)).collect();
        let mut viable_contig_sets = Vec::new();
        for combination in all_combinations {
            let contig_set: Vec<Arc<Contig>> = combination.iter().map(|x| x.bin_contigs.clone().unwrap()).flatten().collect();

            match bin_gen_arc_ref.create_new_bin_from_contigs(&contig_set) {
                Some(result) => {
                    viable_contig_sets.push((contig_set, result.0, result.1, result.2))
                },
                None => continue
            }
            

        }
        viable_contig_sets
    }

    fn filter_eukaryote_and_prokaryote_clusters<'a> (cluster_bin_hashmap: HashMap<usize, Vec<&'a Bin>>) -> HashMap<usize, Vec<&'a Bin>>{
        
        let mut new_cluster_hashmap: HashMap<usize, Vec<&'a Bin>> = HashMap::new();
        let mut cluster_id_count = 0;
        for (cluster, bins) in cluster_bin_hashmap {
            let mut euk_bins = Vec::new();
            let mut prok_bins = Vec::new();
            for bin in bins {
                match bin.bin_type {
                   
                    Some(BinType::eukaryote) => euk_bins.push(bin),
                    Some(BinType::prokaryote) => prok_bins.push(bin),
                    None => panic!("Found bin with no type in clustering - this should not happen.")
                
                }
            }
            if prok_bins.len() > 0 {
        
                new_cluster_hashmap.insert(cluster_id_count as usize, prok_bins);
                cluster_id_count += 1;
       
            }
          
            if euk_bins.len() > 0 {
        
                new_cluster_hashmap.insert(cluster_id_count as usize, euk_bins);
                cluster_id_count += 1;
           
            }
            
        }
        new_cluster_hashmap

    }

    fn create_new_cluster_file_from_eukaryote_and_prokaryote_filtering(output_path: &PathBuf, cluster_bin_hashmap: &HashMap<usize, Vec<&Bin>>) {
        let mut cluster_output_file_string = "bin_id\tcluster_id\n".to_string();

        for (cluster_id, bins) in cluster_bin_hashmap {
            for bin in bins {
            
                let bin_string = &format!("{}\t{}\n", &bin.bin_hash, &cluster_id);
                cluster_output_file_string.push_str(&bin_string);
            
            }
        }
        
        let mut new_cluster_file = File::create(output_path).expect("Can't create eukaryote prokaryote filtered cluster file");

        new_cluster_file.write_all(cluster_output_file_string.as_bytes()).expect("Error could not write cluster string into filtered cluster file");

    }
}   


#[cfg(test)]
mod tests {
    use crate::bin_info_storage::BinType;

    use super::*;
    use std::{fs, env, path::{Path, PathBuf}};
    use lazy_static::lazy_static;

    lazy_static! {
        static ref TEST_CLUSTER_DATA_DIR: PathBuf = PathBuf::from("tests/test_data/cluster_testing/");
    }

    fn create_fake_contigs_for_unit_test() -> Vec<Arc<Contig>> {
        let fake_contig_1 = Contig::new_contig(">fake_contig_1".to_string(), "ATGGCTAGCATCGATGCTAGCAGGAGCGGAGAGCTATGCATGC\n".to_string());
        let fake_contig_2 = Contig::new_contig(">fake_contig_2".to_string(), "ATGGCTAGCATCGATGCTTTTTTAGAGCTATGCATGC\n".to_string());
        let fake_contig_3 = Contig::new_contig(">fake_contig_3".to_string(), "AGTTTTTTGCTAGCATCGATGCTTTTTTGGGGGGC\n".to_string());
        let fake_contig_4 = Contig::new_contig(">fake_contig_4".to_string(), "ATGGCTAGCATCGATGCTTTTTTAGAGCTATGCATGC\n".to_string());
        let fake_contig_5 = Contig::new_contig(">fake_contig_5".to_string(), "ATGGCTAGCAAAATCGATGCTTTTTTAGAGCTATGCATGC\n".to_string());
        let fake_contig_6 = Contig::new_contig(">fake_contig_6".to_string(), "ATGGATATATATCTAGCATCGATGCTTTGTTAGAGCTACCCTGCATGC\n".to_string());
        let fake_contig_7 = Contig::new_contig(">fake_contig_6".to_string(), "ATGGCTACCCGCATCGATGCTTTGTTAGAGCTACCCTGCATGC\n".to_string());
        let fake_contig_8 = Contig::new_contig(">fake_contig_6".to_string(), "ATGGGGGGCTAGCATCGATGCTTTGTTAGAGCTACCCTGCATGC\n".to_string());
        let contig_vec = vec![fake_contig_1, fake_contig_2, fake_contig_3, fake_contig_4, fake_contig_5, fake_contig_6, fake_contig_7, fake_contig_8];
        let fake_contig_arc: Vec<Arc<Contig>> = contig_vec.into_iter().map(|x| Arc::new(x)).collect();
        fake_contig_arc
    }
    fn create_fake_bin_for_unit_test() -> Bin {
       
        
        let mut fake_bin = Bin {fasta_path: PathBuf::from("not_needed"),
        bin_contigs: None,
        completeness: Some(100.0),
        contamination: Some(8.75),
        bin_type: Some(BinType::prokaryote),
        bin_hash: "fake_bin_1".to_string(),
        bin_dir_path: PathBuf::from("not needed")

        };
        fake_bin
    }
    fn create_fake_bin2_for_unit_test() -> Bin {
       
        
        let mut fake_bin2 = Bin {fasta_path: PathBuf::from("not_needed"),
        bin_contigs: None,
        completeness: Some(100.0),
        contamination: Some(9.0),
        bin_type: Some(BinType::prokaryote),
        bin_hash: "fake_bin_2".to_string(),
        bin_dir_path: PathBuf::from("not needed")

        }; 
        fake_bin2
    }
    fn create_fake_bin3_for_unit_test() -> Bin {
        let mut fake_bin3 = Bin {fasta_path: PathBuf::from("not_needed"),
        bin_contigs: None,
        completeness: Some(40.0),
        contamination: Some(9.0),
        bin_type: Some(BinType::prokaryote), 
        bin_hash: "fake_bin_3".to_string(),
        bin_dir_path: PathBuf::from("not needed")

        }; 
        fake_bin3
    }
    fn create_fake_bin4_for_unit_test() -> Bin {
        let mut fake_bin4 = Bin {fasta_path: PathBuf::from("not_needed"),
        bin_contigs: None,
        completeness: Some(40.0),
        contamination: Some(9.0),
        bin_type: Some(BinType::prokaryote),
        bin_hash: "fake_bin_4".to_string(),
        bin_dir_path: PathBuf::from("not needed")

        }; 
        fake_bin4
    }
    fn create_bin_contig_mock() -> (Vec<Bin>, Vec<Arc<Contig>>) {
        let mut bin_1 = create_fake_bin_for_unit_test();
        let mut bin_2 = create_fake_bin2_for_unit_test();
        let mut bin_3 = create_fake_bin3_for_unit_test();
        let mut bin_4 = create_fake_bin4_for_unit_test();

        let fake_contigs = create_fake_contigs_for_unit_test();

        let fake_contigs_for_bin_1 = vec![Arc::clone(&fake_contigs[0]), Arc::clone(&fake_contigs[1]), Arc::clone(&fake_contigs[2]), Arc::clone(&fake_contigs[3]), Arc::clone(&fake_contigs[4])];
        let fake_contigs_for_bin_2 = vec![Arc::clone(&fake_contigs[0]), Arc::clone(&fake_contigs[1]), Arc::clone(&fake_contigs[2]), Arc::clone(&fake_contigs[3])];
        let fake_contigs_for_bin_3 = vec![Arc::clone(&fake_contigs[5]), Arc::clone(&fake_contigs[6]), Arc::clone(&fake_contigs[7])];
        let fake_contigs_for_bin_4 = vec![Arc::clone(&fake_contigs[5]), Arc::clone(&fake_contigs[6]), Arc::clone(&fake_contigs[7])];
        bin_1.bin_contigs = Some(fake_contigs_for_bin_1);
        bin_2.bin_contigs = Some(fake_contigs_for_bin_2);
        bin_3.bin_contigs = Some(fake_contigs_for_bin_3);
        bin_4.bin_contigs = Some(fake_contigs_for_bin_4);
        (vec![bin_1, bin_2, bin_3, bin_4], fake_contigs)
    }


    #[test]
    fn test_contig_presence_vector_gen() {
   
        let fake_contig_arc = create_fake_contigs_for_unit_test();
        let fake_contig_cloned_arcs: Vec<Arc<Contig>> = fake_contig_arc.iter().map(|x| Arc::clone(x)).collect();
        let mut fake_contig_cloned_arcs_2: Vec<Arc<Contig>> = fake_contig_arc.iter().map(|x| Arc::clone(x)).collect();
        fake_contig_cloned_arcs_2.pop();
        let fake_presence_vector = PrepareForClustering::create_contig_presence_vector_for_bin(&fake_contig_cloned_arcs_2, &fake_contig_cloned_arcs);
        assert_eq!(vec![1, 1, 1, 1, 1, 1, 1, 0,], fake_presence_vector);
    
    }
    #[test]
    fn test_create_contig_presence_hashmap() {
        let (mock_bin_vec, fake_contigs) = create_bin_contig_mock();
        let contig_presence_hashmap = PrepareForClustering::create_contig_presence_hashmap(&mock_bin_vec, &fake_contigs);
        let expected_contig_hashmap: HashMap<String, Vec<usize>> = HashMap::from([
            ("fake_bin_1".to_string(), vec![1, 1, 1, 1, 1, 0, 0, 0]),
            ("fake_bin_2".to_string(), vec![1, 1, 1, 1, 0, 0, 0, 0]),
            ("fake_bin_3".to_string(), vec![0, 0, 0, 0, 0, 1, 1, 1]),
            ("fake_bin_4".to_string(), vec![0, 0, 0, 0, 0, 1, 1, 1])
        ]);
        assert_eq!(contig_presence_hashmap, expected_contig_hashmap);

    }
    #[test]
    fn test_presence_file_creation() {
        let mock_contig_hashmap: HashMap<String, Vec<usize>> = HashMap::from([
            ("fake_bin_1".to_string(), vec![1, 1, 1, 1, 1, 0, 0, 0]),
            ("fake_bin_2".to_string(), vec![1, 1, 1, 1, 0, 0, 0, 0]),
            ("fake_bin_3".to_string(), vec![0, 0, 0, 0, 0, 1, 1, 1]),
            ("fake_bin_4".to_string(), vec![0, 0, 0, 0, 0, 1, 1, 1])
        ]);
        let contig_presence_path = &TEST_CLUSTER_DATA_DIR.join("cluster_test_contig_presences.tsv");
        fs::remove_file(&contig_presence_path).unwrap_or(println!("file not present at the moment."));
        PrepareForClustering::write_contig_presence_info_to_file(contig_presence_path, 8, mock_contig_hashmap);
        let mut cluster_gen_file = File::open(&contig_presence_path).unwrap();
        let mut cluster_lines = String::new();
        cluster_gen_file.read_to_string(&mut cluster_lines).unwrap();
        let line_vec: Vec<&str> = cluster_lines.lines().collect();
        if line_vec.len() != 5 {
            panic!("presence file badly created")
        }
        for line in cluster_lines.lines() {
            let tab_count = line.matches("\t").count();
            if tab_count != 8 {
                panic!("presence file badly created")
            }
        }
       
    

    }
    #[test]
    fn test_python_clustering_script() {
       test_presence_file_creation();
       let contig_presence_path = &TEST_CLUSTER_DATA_DIR.join("cluster_test_contig_presences.tsv");
       let cluster_output_path = &TEST_CLUSTER_DATA_DIR.join("cluster_test_output.tsv");
       fs::remove_file(&cluster_output_path).unwrap_or(println!("file not present at the moment."));
       Clustering::run_python_clustering_script(contig_presence_path, cluster_output_path);
       assert!(cluster_output_path.is_file());

    }
    #[test]
    fn test_cluster_bins() {
        let (fake_bins, fake_contigs) = create_bin_contig_mock();
        let contig_presence_path = &TEST_CLUSTER_DATA_DIR.join("cluster_test_contig_presences.tsv");
        let cluster_output_path = &TEST_CLUSTER_DATA_DIR.join("cluster_test_output.tsv");
        let test_hashmap = Clustering::cluster_bins(contig_presence_path, cluster_output_path);
        let expected_hashmap: HashMap<String, usize> = HashMap::from([
            ("fake_bin_1".to_string(), 0),
            ("fake_bin_2".to_string(), 0),
            ("fake_bin_3".to_string(), 1),
            ("fake_bin_4".to_string(), 1),
        ]);
        assert_eq!(test_hashmap, expected_hashmap);
        let result = Clustering::generate_cluster_id_to_clustered_bins_hashmap(test_hashmap, &fake_bins);
        println!("{:?}", result);
        
   //     let expected_bin_hashmap: HashMap<usize, Vec<&Bin>> = HashMap::from([(0, vec![&fake_bins[0], &fake_bins[1]]), (1, vec![&fake_bins[2], &fake_bins[3]])
     //   ]);
       // assert_eq!(expected_bin_hashmap, result);

    }

}