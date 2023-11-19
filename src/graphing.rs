
use itertools::Itertools;
use log::info;
use petgraph::graph::Node;
use petgraph::visit::NodeRef;
use petgraph::{stable_graph::NodeIndex, visit::Bfs};
use rayon::prelude::{IntoParallelIterator, ParallelIterator, IntoParallelRefIterator};
use std::collections::{HashMap, HashSet};
use std::path::PathBuf;
use std::sync::Arc;
use petgraph::{Graph, Undirected, visit::IntoNodeReferences};

use crate::bin_generator::{BinGen, BinGenerator};
use crate::bin_info_storage::{BinInfoStorage, Bin, BinType};
use crate::bin_sets::BinSet;
use crate::contigs::Contig;


pub fn run_graph_clustering(the_bins: Vec<Bin>, bin_generator: Arc<BinGen>, minimum_similarity: f64, cluster_output_directory: PathBuf)  -> Vec<(Vec<Arc<Contig>>, String, f64, f64)>{
    
    let unique_bins = ClusteringPrep::remove_duplicate_bins(the_bins);
    let connected_bins = ClusteringPrep::get_bins_with_minimum_similarity_jaccard_shared(&unique_bins, minimum_similarity);
    let bin_distance_graph = BinDistanceGraph::generate_graph(unique_bins.clone(), connected_bins);
    let new_bin_finder = NewBinFinder{};
    let arc_bin_graph = Arc::new(bin_distance_graph);
    let successful_bins = new_bin_finder.test_each_node(arc_bin_graph, Arc::clone(&bin_generator));
    let all_successful_bins = successful_bins.into_iter().collect_vec();
    let mut all_bins_as_bins: Vec<Bin> = all_successful_bins.into_iter().map(|contig_set| bin_generator.generate_new_bin_from_contigs(contig_set).unwrap()).collect();
    info!("Graph clustering successfully generated {} hybrid bins from {} unique bins", all_bins_as_bins.len(), unique_bins.len());

    let eukaryotic_bins = all_bins_as_bins.iter()
        .filter(|bin| bin.bin_type == BinType::eukaryote).map(|bin| bin.clone()).collect_vec();
    
    let all_eukaryotic_bin_pairs = ClusteringPrep::get_all_eukaryotic_bin_pairs(&eukaryotic_bins);
    let bin_distance_graph = BinDistanceGraph::generate_graph(eukaryotic_bins.clone(), all_eukaryotic_bin_pairs);
    let new_bin_finder = NewBinFinder{};
    let arc_bin_graph = Arc::new(bin_distance_graph);
    let euk_successful_bins = new_bin_finder.test_each_node(arc_bin_graph, Arc::clone(&bin_generator));
    let all_euk_successful_bins: Vec<Bin> = euk_successful_bins.into_iter().map(|contig_set| bin_generator.generate_new_bin_from_contigs(contig_set).unwrap()).collect();
    all_bins_as_bins.extend(all_euk_successful_bins);

    let all_extended_unique_bins = ClusteringPrep::remove_duplicate_bins(all_bins_as_bins).into_iter().map(|bin| Arc::new(bin)).collect_vec();
    let bin_set_of_bins_produced_by_clustering = BinSet::make_bin_set_from_bins_vec(all_extended_unique_bins);
    bin_set_of_bins_produced_by_clustering.create_bin_set_dir_and_info_from_best_hashes(&bin_generator.hash_directory, &cluster_output_directory, false);
    let all_bins_by_information = bin_set_of_bins_produced_by_clustering.bins.iter()
        .map(|bin| (bin.bin_contigs.clone(), bin.bin_hash.clone(), bin.completeness.clone(), bin.contamination.clone())).collect_vec();

   // let connected_routes = new_bin_finder.find_connected_routes(&bin_distance_graph);
    // let all_successful_bins = new_bin_finder.test_each_connected_route(bin_distance_graph, connected_routes, bin_generator);
    info!("Graph clustering additional eukaryotic stage: successfully generated {} hybrid bins from {} unique bins", all_bins_by_information.len(), unique_bins.len());
    info!("Graph clustering stage complete!");
    all_bins_by_information
 
}
struct ClusteringPrep;


impl ClusteringPrep {
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
/* 
    fn get_all_bin_pairs<'a>(bins: &'a Vec<Bin>, minimum_similarity: f64) -> Vec<(&'a Bin, &'a Bin)> {
        let eukaryotic_bins = bins.iter()
            .filter(|bin| bin.bin_type == BinType::eukaryote).collect_vec();
        
        let eukaryotic_bin_combinations = ClusteringPrep::get_all_eukaryotic_bin_pairs(bins);
        let prokaryotic_bins = bins.iter().filter(|bin| bin.bin_type == BinType::prokaryote).collect_vec();
        let prokaryotic_bin_combinations = ClusteringPrep::get_bins_with_minimum_similarity_jaccard_shared(bins, minimum_similarity)


    } */
    fn get_all_eukaryotic_bin_pairs<'a>(bins: &'a Vec<Bin>) -> Vec<(&'a Bin, &'a Bin)> {

        let all_potential_eukaryotic_bin_pairs: Vec<(&'a Bin, &'a Bin)> = bins.into_iter()
            .combinations(2)
            .filter_map(|bin_combination| ClusteringPrep::check_whether_pair_might_improve(bin_combination[0], bin_combination[1]))
            .collect();
        
        all_potential_eukaryotic_bin_pairs
    }
    fn check_whether_pair_might_improve<'a>(bin_1: &'a Bin, bin_2: &'a Bin) -> Option<(&'a Bin, &'a Bin)> { 
        let bin_1_contigs = bin_1.bin_contigs.iter().collect_vec();
        let bin_2_contigs = bin_2.bin_contigs.iter().collect_vec();
        
        let mut bin_1_unique_contig_markers: HashSet<String> = bin_1_contigs.iter() // markers that are in contigs that only bin 1 has
            .filter(|contig| !bin_2_contigs.contains(contig))
            .filter_map(|contig| contig.eukaryotic_contig_info.clone())
            .map(|contig_euk_info| contig_euk_info.complete_buscos)
            .flatten().collect();

        let mut bin_2_unique_contig_markers: HashSet<String> = bin_2_contigs.iter() // markers that are in contigs that only bin 1 has
            .filter(|contig| !bin_1_contigs.contains(contig))
            .filter_map(|contig| contig.eukaryotic_contig_info.clone())
            .map(|contig_euk_info| contig_euk_info.complete_buscos)
            .flatten().collect();
        // logic here is it checks if combining the contigs unique to each bin will result in an increase in unique markers (unique meaning increased completeness)
        let combined_bin_1_2_unique_markers: HashSet<&String> = bin_1_unique_contig_markers.iter().chain(bin_2_unique_contig_markers.iter()).collect();
        if (combined_bin_1_2_unique_markers.len() > bin_1_unique_contig_markers.len()) & (combined_bin_1_2_unique_markers.len() > bin_2_unique_contig_markers.len()) {
            Some((bin_1, bin_2))
        } else {
            None
        }


        // let bin_1_eukaryotic_markers = bin_1.bin_contigs.iter()
        //    .filter_map(|contig| contig.eukaryotic_contig_info.clone())
         //   .map(|euk_contig_info| euk_contig_info.complete_buscos)
         //   .collect_vec();
    }


    fn get_bins_with_minimum_similarity_jaccard_shared<'a>(bins: &'a Vec<Bin>, minimum_similarity: f64) -> Vec<(&'a Bin, &'a Bin)>{
        let mut bins_with_minimum_similarity = Vec::new();
        let all_potential_bin_pairs: Vec<Vec<&'a Bin>> = bins.into_iter().combinations(2).collect();
        for bin_pair in all_potential_bin_pairs {

            if ClusteringPrep::calc_jaccard_similarity(bin_pair[0], bin_pair[1]) > minimum_similarity {

                    bins_with_minimum_similarity.push((bin_pair[0], bin_pair[1]));

                }

        }
        bins_with_minimum_similarity
    }

    fn calc_jaccard_similarity(bin_1: &Bin, bin_2: &Bin) -> f64 {
        let bin_1_contigs = bin_1.bin_contigs.clone();
        let bin_2_contigs = bin_2.bin_contigs.clone();

        let bin_1_contigs_set: HashSet<_> = bin_1_contigs.iter().collect();
        let bin_2_contigs_set: HashSet<_> = bin_2_contigs.iter().collect();

        let intersection_size = bin_1_contigs_set.intersection(&bin_2_contigs_set).count();
        let union_size = bin_1_contigs_set.union(&bin_2_contigs_set).count();

        if union_size == 0 {
            return 0.0;
        }

        intersection_size as f64 / union_size as f64
    }


}



#[derive(Clone)]
pub struct BinDistanceGraph  {
    the_graph: Graph::<Bin, f64>,
    node_bin_dict: HashMap<NodeIndex, Bin>,

}

impl BinDistanceGraph {

    fn generate_graph(list_of_initial_bins: Vec<Bin>, list_of_connected_bins: Vec<(&Bin, &Bin)>) -> BinDistanceGraph {
        let mut graph = Graph::<Bin, f64>::new();
        let node_bin_dict = BinDistanceGraph::add_nodes_to_graph(&mut graph, list_of_initial_bins);
        let mut the_graph = BinDistanceGraph {
            the_graph: graph,
            node_bin_dict: node_bin_dict
        };
        the_graph.add_edges_to_graph(list_of_connected_bins);
        the_graph
    }

    fn add_nodes_to_graph(graph: &mut Graph::<Bin, f64>, list_of_initial_bins: Vec<Bin> ) -> HashMap<NodeIndex, Bin> {

        let mut bin_nodes = Vec::new();
        let mut node_bin_dict = HashMap::new();
        println!("bins for add node to graph: {}", list_of_initial_bins.len());
        for bin in list_of_initial_bins.into_iter() {
            let node = graph.add_node(bin.clone());
            bin_nodes.push(node);

            node_bin_dict.insert(node, bin);

        }

        node_bin_dict

    }

    fn add_edges_to_graph(&mut self, connected_bins: Vec<(&Bin, &Bin)>) {


        for (bin_1, bin_2) in connected_bins {
            let node_1 = self.the_graph.node_indices().find(|&i| &self.the_graph[i] == bin_1).unwrap();
            let node_2 = self.the_graph.node_indices().find(|&i| &self.the_graph[i] == bin_2).unwrap();
            self.the_graph.add_edge(node_1, node_2, 1.0);

        }

    }


}

struct NewBinFinder;

impl NewBinFinder {


    fn test_each_node(&self, bin_distance_graph: Arc<BinDistanceGraph>, bin_gen_arc: Arc<BinGen>, ) -> HashSet<Vec<Arc<Contig>>> {

        let unique_created_bins: HashSet<Vec<Arc<Contig>>> = bin_distance_graph.node_bin_dict.clone().par_iter().map(|(node, bin)| {
            let mut new_bins = vec![bin.bin_contigs.clone()];
            self.test_node_potential_bins(&Arc::clone(&bin_distance_graph), vec![node], &Arc::clone(&bin_gen_arc), (bin.completeness, bin.contamination), &mut new_bins);
            info!("One node finished!");
            new_bins
        }).flatten().collect();
      //  let mut unique_bins = HashSet::new();
        unique_created_bins
    }


    fn test_node_potential_bins(&self, bin_distance_graph: &BinDistanceGraph, current_bin_nodes: Vec<&NodeIndex>, bin_generator: &Arc<BinGen>, current_bin_quality: (f64, f64), successful_bins: &mut Vec<Vec<Arc<Contig>>>) {
        let current_bin_test_node = current_bin_nodes[current_bin_nodes.len() - 1];
        let neighbor_nodes: Vec<NodeIndex> = bin_distance_graph.the_graph.neighbors_undirected(current_bin_test_node.clone())
            .filter(|neighbor_node| !current_bin_nodes.contains(&neighbor_node)).collect();

        let current_bins: Vec<&Bin> = current_bin_nodes.iter().map(|node| bin_distance_graph.node_bin_dict.get(&node).unwrap()).collect();
        
        for neighbor_node in neighbor_nodes {

            let hypothetical_bin = bin_distance_graph.node_bin_dict.get(&neighbor_node).unwrap();
            let mut bin_test = current_bins.clone();
            bin_test.push(hypothetical_bin);
            let contigs_to_test: Vec<Arc<Contig>> = bin_test.iter().map(|bin| bin.bin_contigs.clone()).flatten().unique().collect();
            let mut single_vec = Vec::new();
            let mut intersect_contig_test = Vec::new();

            bin_test.iter().map(|bin| bin.bin_contigs.clone())
                .for_each(|x| {
                    if !single_vec.contains(&x) {
                        intersect_contig_test.extend(x);
                    } else {
                        single_vec.push(x)
                    }
                });

            match bin_generator.generate_new_bin_from_contigs(intersect_contig_test.clone()) {
                Some(bin_res) => {
                    if bin_res.completeness > current_bin_quality.0 || bin_res.contamination < current_bin_quality.1 {
                        let mut current_bin_nodes_plus_successful_neighbor = current_bin_nodes.clone();
                        current_bin_nodes_plus_successful_neighbor.push(&neighbor_node);
                        successful_bins.push(intersect_contig_test);
                        self.test_node_potential_bins(bin_distance_graph, current_bin_nodes_plus_successful_neighbor,  &bin_generator, (bin_res.completeness, bin_res.contamination), successful_bins);
                    }
                },
                None => ()
            }

        
            

            match bin_generator.generate_new_bin_from_contigs(contigs_to_test.clone()) {
                Some(bin_res) => {
                    if bin_res.completeness > current_bin_quality.0 || bin_res.contamination < current_bin_quality.1 {

                        let mut current_bin_nodes_plus_successful_neighbor = current_bin_nodes.clone();
                        current_bin_nodes_plus_successful_neighbor.push(&neighbor_node);
                        successful_bins.push(contigs_to_test);
                        self.test_node_potential_bins(bin_distance_graph, current_bin_nodes_plus_successful_neighbor,  &bin_generator, (bin_res.completeness, bin_res.contamination), successful_bins);
                    }
                },
                None => ()
            }

        }

    }



}


/* 
mod tests {

    use crate::{bin_info_storing::BinSet, bin_classes::BinType, initial_bins_and_contigs::gather_initial_bins_and_contig_information};

    use super::*;
    use std::{fs, env, path::{Path, PathBuf}};
    use lazy_static::lazy_static;
    lazy_static! {
    static ref TEST_DATA_HASH: PathBuf = PathBuf::from("tests/test_data/graphing_testing/hash_directory/");
}
    lazy_static! {
    static ref TEST_DATA_DIR: PathBuf = PathBuf::from("tests/test_data/graphing_testing/test_bins/");
    }
    lazy_static! { // unit test will only work if database is added here, not uploaded to git repo due to size
        static ref COMPLEASM_DB_LIB: PathBuf = PathBuf::from("tests/test_data/databases_for_testing/");
    }
    lazy_static! { // unit test will only work if database is added here, not uploaded to git repo due to size

        static ref CHECKM2_DB_PATH: PathBuf = PathBuf::from("tests/test_data/databases_for_testing/uniref100.KO.1.dmnd");

    }




    fn create_fake_contigs_for_unit_test() -> Vec<Arc<Contig>> {
        let fake_contig_1 = Contig::new_contig(">fake_contig_1".to_string(), "ATGGCTAGCATCGATGCTAGCAGGAGCGGAGAGCTATGCATGC\n".to_string());
        let fake_contig_2 = Contig::new_contig(">fake_contig_2".to_string(), "ATGGCTAGCATCGATGCTTTTTTAGAGCTATGCATGC\n".to_string());
        let fake_contig_3 = Contig::new_contig(">fake_contig_3".to_string(), "AGTTTTTTGCTAGCATCGATGCTTTTTTGGGGGGC\n".to_string());
        let fake_contig_4 = Contig::new_contig(">fake_contig_4".to_string(), "ATGGCTAGCATCGATGCTTTTTTAGAGCTATGCATGC\n".to_string());
        let fake_contig_5 = Contig::new_contig(">fake_contig_5".to_string(), "ATGGCTAGCAAAATCGATGCTTTTTTAGAGCTATGCATGC\n".to_string());
        let fake_contig_6 = Contig::new_contig(">fake_contig_6".to_string(), "ATGGATATATATCTAGCATCGATGCTTTGTTAGAGCTACCCTGCATGC\n".to_string());
        let fake_contig_7 = Contig::new_contig(">fake_contig_7".to_string(), "ATGGCTACCCGCATCGATGCTTTGTTAGAGCTACCCTGCATGC\n".to_string());
        let fake_contig_8 = Contig::new_contig(">fake_contig_8".to_string(), "ATGGGGGGCTAGCATCGATGCTTTGTTAGAGCTACCCTGCATGC\n".to_string());
        let contig_vec = vec![fake_contig_1, fake_contig_2, fake_contig_3, fake_contig_4, fake_contig_5, fake_contig_6, fake_contig_7, fake_contig_8];
        let fake_contig_arc: Vec<Arc<Contig>> = contig_vec.into_iter().map(|x| Arc::new(x)).collect();
        fake_contig_arc
    }
    fn create_fake_contig_set_bins() -> Vec<(Vec<Arc<Contig>>, String, f64, f64)> {
        let fake_contigs = create_fake_contigs_for_unit_test();
        let fake_first_best_contigs = vec![Arc::clone(&fake_contigs[0]), Arc::clone(&fake_contigs[1]), Arc::clone(&fake_contigs[2]), Arc::clone(&fake_contigs[7])]; // has one the same as last, zero the same as second
        let fake_first_best = (fake_first_best_contigs, "fake_hsh_1".to_string(), 99.9, 0.00);
        let fake_second_best_contigs = vec![Arc::clone(&fake_contigs[3]), Arc::clone(&fake_contigs[4]), Arc::clone(&fake_contigs[5]), Arc::clone(&fake_contigs[6])];
        let fake_second_best = (fake_second_best_contigs, "fake_hsh_2".to_string(), 90.0, 1.0);
        let fake_worst_contigs = vec![Arc::clone(&fake_contigs[0]), Arc::clone(&fake_contigs[3]), Arc::clone(&fake_contigs[5]), Arc::clone(&fake_contigs[6]), Arc::clone(&fake_contigs[4])]; // has one contig the same as the first, four the same as the second
        let fake_worst = (fake_worst_contigs, "fake_hsh_3".to_string(), 10.0, 100.0);

        vec![fake_first_best, fake_second_best, fake_worst]

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

        let fake_contigs = create_fake_contigs_for_unit_test();

        let fake_contigs_for_bin_1 = vec![Arc::clone(&fake_contigs[0]), Arc::clone(&fake_contigs[1]), Arc::clone(&fake_contigs[2]), Arc::clone(&fake_contigs[7])]; // has one the same as last, zero the same as second
        let fake_contigs_for_bin_2 = vec![Arc::clone(&fake_contigs[3]), Arc::clone(&fake_contigs[4]), Arc::clone(&fake_contigs[5]), Arc::clone(&fake_contigs[6])];
        let fake_contigs_for_bin_3 = vec![Arc::clone(&fake_contigs[0]), Arc::clone(&fake_contigs[3]), Arc::clone(&fake_contigs[5]), Arc::clone(&fake_contigs[6]), Arc::clone(&fake_contigs[4])]; // has one contig the same as the first, four the same as the second
        bin_1.bin_contigs = Some(fake_contigs_for_bin_1);
        bin_2.bin_contigs = Some(fake_contigs_for_bin_2);
        bin_3.bin_contigs = Some(fake_contigs_for_bin_3);
        (vec![bin_1, bin_2, bin_3], fake_contigs)
    }


        #[test]
        fn test_jaccard_similarity() {
            let (fake_bins, _) = create_bin_contig_mock();
            let bin_one_against_two_result = 0.0;
            let bin_one_against_three_result = 0.125;
            let bin_two_against_three_result = 0.8;
            assert_eq!(ClusteringPrep::calc_jaccard_similarity(&fake_bins[0], &fake_bins[1]), bin_one_against_two_result);
            assert_eq!(ClusteringPrep::calc_jaccard_similarity(&fake_bins[0], &fake_bins[2]), bin_one_against_three_result);
            assert_eq!(ClusteringPrep::calc_jaccard_similarity(&fake_bins[1], &fake_bins[2]), bin_two_against_three_result);
            let connected_bins = ClusteringPrep::get_bins_with_minimum_similarity_jaccard_shared(&fake_bins, 0.7);
            assert_eq!(connected_bins.len(), 1);
            assert!(connected_bins[0] == (&fake_bins[1], &fake_bins[2]));

        }
        #[test]
        fn test_graph_creation() {
            let (fake_bins, _) = create_bin_contig_mock();
            let connected_bins = ClusteringPrep::get_bins_with_minimum_similarity_jaccard_shared(&fake_bins, 0.7);
            let bin_distance_graph = BinDistanceGraph::generate_graph(fake_bins.clone(), connected_bins);
            let test_graph = bin_distance_graph.the_graph;
            let mut fake_graph = Graph::<Bin, f64>::new();

            let f_node0 = fake_graph.add_node(fake_bins[0].clone());
            let f_node1 = fake_graph.add_node(fake_bins[1].clone());
            let f_node2 = fake_graph.add_node(fake_bins[2].clone());

            fake_graph.add_edge(f_node1, f_node2, 1.0);
            if fake_graph.node_count() != test_graph.node_count() {
                panic!("Node count different.")
            }

            // Check edge count
            if fake_graph.edge_count() != test_graph.edge_count() {
                panic!("Edge count different.")
            }

            // Check nodes and node weights
            for node in fake_graph.node_indices() {
                if fake_graph[node] != test_graph[node] {
                    panic!("Different nodes")
                }
            }


            for edge in fake_graph.edge_indices() {
                let (a, b) = fake_graph.edge_endpoints(edge).unwrap();
                let weight1 = &fake_graph[edge];
                match test_graph.find_edge(a, b) {
                    Some(edge2) => {
                        if weight1 != &test_graph[edge2] {
                            panic!("Edge weight does not match")
                        }
                    }
                    None => panic!("Cannot find edge")
                }
            }

        }
        #[test]
        fn real_example_bins(){
            
            let bin_generator = BinGen::initialise_bin_gen(Some(CHECKM2_DB_PATH.to_path_buf()), Some(COMPLEASM_DB_LIB.to_path_buf()), TEST_DATA_HASH.clone(), 100.0, 0.0);
            let arc_bin_gen = Arc::new(bin_generator);
            let (bins, contigs) = gather_initial_bins_and_contig_information(&TEST_DATA_DIR, Arc::clone(&arc_bin_gen));
            let connected_bins = ClusteringPrep::get_bins_with_minimum_similarity_jaccard_shared(&bins, 0.1);
            println!("Number of connected bin pairs equals: {}", connected_bins.len());
            for (bin0, bin1) in &connected_bins {
                if bin0.bin_hash == "08c929dd57d754fe8e56f2befc8db4daed8bded44d54525e970e9be228345b9c" || bin1.bin_hash == "08c929dd57d754fe8e56f2befc8db4daed8bded44d54525e970e9be228345b9c" {
                    panic!("Found bin that should not be connected to others.")
                }
            }
            
            let bin_dist_graph = BinDistanceGraph::generate_graph(bins.clone(), connected_bins);
            let new_bin_finder = NewBinFinder{};
            let test_node_hash = "e18f3f062e2f065e9792922345b0948e5449d627e6a3c5e302d69939a53109e6".to_string();
            let mut test_node = None;
            let mut test_node_contigs = None;
            for (node, bin) in &bin_dist_graph.node_bin_dict {
                println!("hash: {}", bin.bin_hash);
                if bin.bin_hash == test_node_hash {
                    println!("found");
                    test_node = Some(node);
                    test_node_contigs = Some(bin.bin_contigs.clone());
                    break
                }
            }

      //      bin_dist_graph.
            let mut generated_bins = vec![test_node_contigs.unwrap().unwrap()];
            new_bin_finder.test_node_potential_bins(&bin_dist_graph, vec![&test_node.unwrap()], &Arc::clone(&arc_bin_gen), (91.42, 9.81), &mut generated_bins);
            assert_eq!(generated_bins.len(), 3);
            let all_bins = new_bin_finder.test_each_node(Arc::new(bin_dist_graph), Arc::clone(&arc_bin_gen));
            println!("{}", all_bins.len());
            panic!("done");
            let test_connected_routes = new_bin_finder.find_connected_routes(&bin_dist_graph);
            assert_eq!(test_connected_routes.len(), 2);

           
            if test_connected_routes[0].len() == 1 || test_connected_routes[1].len() == 1 && (test_connected_routes[0].len() != test_connected_routes[1].len()) {
           
                panic!("Find connected routes method failure: Couldn't find the route with one node");
           
            }
           
            for route in &test_connected_routes {
           
                if route.len() == 1 {
           
                    let bin = bin_dist_graph.node_bin_dict.get(&route[0]).clone().unwrap();
           
                    if !(bin.bin_hash == "08c929dd57d754fe8e56f2befc8db4daed8bded44d54525e970e9be228345b9c") {
           
                        panic!("find connected routes method failure: wrong bin hash in route with one node");
                    }
                }
            }
            let total_potential_bins: Vec<Vec<Arc<Contig>>> = new_bin_finder.test_each_connected_route(bin_dist_graph, test_connected_routes, Arc::clone(&arc_bin_gen));;
           
            println!("{}", total_potential_bins.len());
          //  new_bin_finder.test_connected_route(nodes_in_route, Arc::new(bin_dist_graph), Arc::clone(&arc_bin_gen));

        }

        #[test]
        fn validate_test_node_potential_bins_method() {
            let new_bin_finder = NewBinFinder{};
          //  new_bin_finder.test_node_potential_bins(bin_distance_graph, current_bin_nodes, bin_generator, current_bin_quality, successful_bins)
        }


}
*/