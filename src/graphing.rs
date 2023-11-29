
use itertools::Itertools;
use log::{info, debug};
use petgraph::graph::Node;
use petgraph::visit::NodeRef;
use petgraph::{stable_graph::NodeIndex, visit::Bfs};
use rayon::iter::IndexedParallelIterator;
use rayon::prelude::{IntoParallelIterator, ParallelIterator, IntoParallelRefIterator};
use std::collections::{HashMap, HashSet};
use std::path::PathBuf;
use std::sync::Arc;
use petgraph::{Graph, Undirected, visit::IntoNodeReferences};

use crate::bin_generator::{BinGen, BinGenerator};
use crate::bin_info_storage::{BinInfoStorage, Bin, BinType};
use crate::bin_sets::BinSet;
use crate::contigs::Contig;


pub fn run_graph_clustering(the_bins: Vec<Bin>, bin_generator: Arc<BinGen>, max_jaccard_distance: f64, max_euclidean_distance: f64, cluster_output_directory: PathBuf)  -> BinSet {
    
    let initial_unique_bins = ClusteringPrep::remove_duplicate_bins(the_bins);
    let connected_bins = ClusteringPrep::get_bin_pairs_with_less_than_max_jaccard_distance(&initial_unique_bins, max_jaccard_distance);
    let bin_distance_graph = BinDistanceGraph::generate_graph(initial_unique_bins.clone(), connected_bins);
    let new_bin_finder = NewBinFinder{};
    let arc_bin_graph = Arc::new(bin_distance_graph);
    let successful_bins = new_bin_finder.test_each_node(arc_bin_graph, Arc::clone(&bin_generator));
    let all_successful_bins = successful_bins.into_iter().collect_vec();
    let mut all_bins_as_bins: Vec<Bin> = all_successful_bins.into_iter().map(|contig_set| bin_generator.generate_new_bin_from_contigs(contig_set).unwrap()).collect();
    info!("Graph clustering successfully generated {} hybrid bins from {} unique bins", all_bins_as_bins.len(), initial_unique_bins.len());



    let all_extended_unique_bins = ClusteringPrep::remove_duplicate_bins(all_bins_as_bins).into_iter().map(|bin| Arc::new(bin)).collect_vec();
    let bin_set_of_bins_produced_by_clustering = BinSet::make_bin_set_from_bins_vec(all_extended_unique_bins);
    bin_set_of_bins_produced_by_clustering.create_bin_set_dir_and_info_from_best_hashes(&bin_generator.hash_directory, &cluster_output_directory, false);
    bin_set_of_bins_produced_by_clustering
 
}


pub fn run_additional_eukaryotic_clustering_stage(bins: &Vec<Bin>, bin_generator: Arc<BinGen>, max_euclidean_distance: f64, cluster_output_directory: PathBuf) -> BinSet{
    
    let eukaryotic_bins = bins.iter()
        .filter(|bin| bin.bin_type == BinType::eukaryote)
        .map(|bin| bin.clone())
        .collect_vec();

    let all_eukaryotic_bin_pairs = ClusteringPrep::get_euk_bin_pairs_with_less_than_max_euclidean_distance(&eukaryotic_bins, 4, max_euclidean_distance);
    info!("Identified {} pairs with minimum_distance", all_eukaryotic_bin_pairs.len());
    let bin_distance_graph = BinDistanceGraph::generate_graph(eukaryotic_bins.clone(), all_eukaryotic_bin_pairs);
    let new_bin_finder = NewBinFinder{};
    let arc_bin_graph = Arc::new(bin_distance_graph);
    let euk_successful_bins = new_bin_finder.test_each_node(arc_bin_graph, Arc::clone(&bin_generator));
    let all_euk_successful_bins: Vec<Bin> = euk_successful_bins.into_iter().map(|contig_set| bin_generator.generate_new_bin_from_contigs(contig_set).unwrap()).collect();
    let unique_euk_bins = ClusteringPrep::remove_duplicate_bins(all_euk_successful_bins).into_iter().map(|bin| Arc::new(bin)).collect_vec();
    let bin_set_of_bins_produced_by_clustering = BinSet::make_bin_set_from_bins_vec(unique_euk_bins);

    bin_set_of_bins_produced_by_clustering.create_bin_set_dir_and_info_from_best_hashes(&bin_generator.hash_directory, &cluster_output_directory, false);
    info!("Graph clustering additional eukaryotic stage: successfully generated {} hybrid bins from {} eukaryotic bins", bin_set_of_bins_produced_by_clustering.bins.len(), eukaryotic_bins.len());
    info!("Graph clustering stage complete!");
    bin_set_of_bins_produced_by_clustering

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

    fn get_bin_pairs_with_less_than_max_jaccard_distance<'a>(bins: &'a Vec<Bin>, max_distance: f64) -> Vec<(&'a Bin, &'a Bin)>{
        let mut bins_with_minimum_similarity = Vec::new();
        let all_potential_bin_pairs: Vec<Vec<&'a Bin>> = bins.into_iter().combinations(2).collect();
        for bin_pair in all_potential_bin_pairs {

            if ClusteringPrep::calc_jaccard_distance(bin_pair[0], bin_pair[1]) < max_distance {

                    bins_with_minimum_similarity.push((bin_pair[0], bin_pair[1]));

                }

        }
        bins_with_minimum_similarity
    }

    fn calc_jaccard_distance(bin_1: &Bin, bin_2: &Bin) -> f64 {
        let bin_1_contigs = bin_1.bin_contigs.clone();
        let bin_2_contigs = bin_2.bin_contigs.clone();

        let bin_1_contigs_set: HashSet<_> = bin_1_contigs.iter().collect();
        let bin_2_contigs_set: HashSet<_> = bin_2_contigs.iter().collect();

        let intersection_size = bin_1_contigs_set.intersection(&bin_2_contigs_set).count();
        let union_size = bin_1_contigs_set.union(&bin_2_contigs_set).count();

        if union_size == 0 {
            return 0.0;
        }

        1.0 - (intersection_size as f64 / union_size as f64)
    }

    fn get_euk_bin_pairs_with_less_than_max_euclidean_distance<'a>(bins: &'a [Bin], kmer_size: usize, max_euclidean_distance: f64) -> Vec<(&'a Bin, &'a Bin)>{
        let mut all_potential_bin_pairs: Vec<(&Bin, &Bin)> = bins.into_iter()
            .combinations(2).map(|bin_combo| (bin_combo[0], bin_combo[1])).collect();
        let all_unique_contigs = bins.iter().map(|bin| bin.bin_contigs.clone()).flatten().unique().collect_vec();
        let contig_kmer_hashmap = ClusteringPrep::create_contig_kmer_dict_from_bins(kmer_size, all_unique_contigs);
        let bin_kmer_dict = ClusteringPrep::create_bin_kmer_dict_from_contig_kmer_hashmap(bins, contig_kmer_hashmap);
        let all_euclidean_pairing_result_sorted: Vec<((&&Bin, &&Bin), f64)> = all_potential_bin_pairs
            .par_iter()
            .map(|(bin_1, bin_2)| ((bin_1, bin_2), ClusteringPrep::calculate_euclidean_distance_between_bin_pair(bin_kmer_dict.get(bin_1).unwrap(), bin_kmer_dict.get(bin_2).unwrap())))
            .collect::<Vec<_>>() // Collect into a vector first
            .into_iter()
            .sorted_by(|a, b| a.1.partial_cmp(&b.1).unwrap()) // Sort the vector by the second item of the tuple
            .collect(); // Finally, collect into a HashMap
        let max_dist = all_euclidean_pairing_result_sorted[0].1;
        let min_dist = all_euclidean_pairing_result_sorted[all_euclidean_pairing_result_sorted.len()].1;

        let mut viable_bin_pairs = Vec::new();
        for (bin_pair, euclid_dist) in all_euclidean_pairing_result_sorted {
            let min_max_normalised_distance = (euclid_dist - min_dist) / (max_dist - min_dist);
            debug!("normalised euclidean distance equals: {}", min_max_normalised_distance);
            if min_max_normalised_distance < max_euclidean_distance {
                viable_bin_pairs.push((*bin_pair.0, *bin_pair.1));
            }
        }
        viable_bin_pairs

    }

    fn create_contig_kmer_dict_from_bins(kmer_size: usize, contigs: Vec<Arc<Contig>>) -> HashMap<Arc<Contig>, Vec<String>> {
        let mut contig_kmer_dict = HashMap::new();
        for contig in contigs {
            let kmers = contig.get_kmers(kmer_size);
            contig_kmer_dict.insert(contig, kmers);
        }
        contig_kmer_dict
    }

    fn create_bin_kmer_dict_from_contig_kmer_hashmap<'a> (bins: &'a [Bin], contig_kmer_hashmap: HashMap<Arc<Contig>, Vec<String>>) -> HashMap<&Bin, HashMap<String, f64>> {
        let bin_kmer_hashmap: HashMap<&Bin, HashMap<String, f64>> = bins.par_iter()
            .map(|bin| (bin, bin.get_bin_kmers_from_contig_kmer_hashmap(&contig_kmer_hashmap).clone()))
            .map(|(bin, bin_kmer_hashmap)| (bin, ClusteringPrep::convert_bin_kmer_dict_to_bin_kmer_frequency_ratio_dict(bin_kmer_hashmap)))
            .collect();
        debug!("Successfully generated bin kmer hashmap!");
        bin_kmer_hashmap
    }
    fn convert_bin_kmer_dict_to_bin_kmer_frequency_ratio_dict(bin_kmer_dict: HashMap<String, i32>) -> HashMap<String, f64> {
        let mut bin_kmer_frequency_ratio_dict = HashMap::new();
        let total_kmers = bin_kmer_dict.iter().fold(0, |acc, (key, value)| acc + value);
        for (kmer, count) in bin_kmer_dict {
            bin_kmer_frequency_ratio_dict.insert(kmer, (count / total_kmers) as f64);
        }
        bin_kmer_frequency_ratio_dict
    }


    fn calculate_euclidean_distance_between_bin_pair(bin_1_kmer_hash: &HashMap<String, f64>, bin_2_kmer_hash: &HashMap<String, f64>) -> f64 {
        // uses tetranucleotide frequency ratio (as a means to normalise contig sizes)
        let all_unique_kmers: HashSet<&String> = bin_1_kmer_hash.keys().chain(bin_2_kmer_hash.keys()).unique().collect();

        let mut sum_of_squared_differences = 0.0;

        for kmer in all_unique_kmers {
            let bin_1_num = *bin_1_kmer_hash.get(kmer).unwrap_or(&0.0);
            let bin_2_num = *bin_2_kmer_hash.get(kmer).unwrap_or(&0.0);
    
            let squared_kmer_diff = (bin_1_num - bin_2_num) * (bin_1_num - bin_2_num);
            sum_of_squared_differences += squared_kmer_diff;
        }
        let euclidean_distance = sum_of_squared_differences.sqrt();
        debug!("Calculated euclidean distance as: {}", euclidean_distance);
        euclidean_distance
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

            let current_bin_contigs: HashSet<Arc<Contig>> = current_bins.clone().into_iter().map(|bin| bin.bin_contigs.clone()).flatten().unique().collect();
            let hypothetical_bin_contigs: HashSet<Arc<Contig>> = hypothetical_bin.bin_contigs.clone().into_iter().collect();
            
            let union_contigs: Vec<Arc<Contig>> = current_bin_contigs.union(&hypothetical_bin_contigs).cloned().collect();
            let intersection_contigs: Vec<Arc<Contig>> = current_bin_contigs.intersection(&hypothetical_bin_contigs).cloned().collect();


            match bin_generator.generate_new_bin_from_contigs(intersection_contigs.clone()) {
                Some(bin_res) => {
                    if self.check_if_improvement_conditions_met(&current_bin_quality, &(bin_res.completeness, bin_res.contamination)) {
                        let mut current_bin_nodes_plus_successful_neighbor = current_bin_nodes.clone();
                        current_bin_nodes_plus_successful_neighbor.push(&neighbor_node);
                        successful_bins.push(intersection_contigs);
                        self.test_node_potential_bins(bin_distance_graph, current_bin_nodes_plus_successful_neighbor,  &bin_generator, (bin_res.completeness, bin_res.contamination), successful_bins);
                    }
                },
                None => ()
            }

        
            

            match bin_generator.generate_new_bin_from_contigs(union_contigs.clone()) {
                Some(bin_res) => {
                    if self.check_if_improvement_conditions_met(&current_bin_quality, &(bin_res.completeness, bin_res.contamination)) {

                        let mut current_bin_nodes_plus_successful_neighbor = current_bin_nodes.clone();
                        current_bin_nodes_plus_successful_neighbor.push(&neighbor_node);
                        successful_bins.push(union_contigs);
                        self.test_node_potential_bins(bin_distance_graph, current_bin_nodes_plus_successful_neighbor,  &bin_generator, (bin_res.completeness, bin_res.contamination), successful_bins);
                    }
                },
                None => ()
            }

        }

    }
    fn check_if_improvement_conditions_met(&self, current_bin_quality: &(f64, f64), new_potential_bin_quality: &(f64, f64)) -> bool {
        
        let change_in_completeness = new_potential_bin_quality.0 - current_bin_quality.0;
        let change_in_contamination = new_potential_bin_quality.0 - current_bin_quality.1;
        change_in_completeness >= change_in_contamination

        
    }



}



mod tests {

    use crate::{bin_info_storage::BinType};

    use super::*;
    use std::{fs, env, path::{Path, PathBuf}};
    use lazy_static::lazy_static;
    use rand::rngs::mock;
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
        let fake_contig_1 = Contig::new_contig(">fake_contig_1".to_string(), "AAAAAA".to_string());
        let fake_contig_2 = Contig::new_contig(">fake_contig_2".to_string(), "TTTTTT".to_string());
        let fake_contig_3 = Contig::new_contig(">fake_contig_3".to_string(), "CCCCCC".to_string());
        let fake_contig_4 = Contig::new_contig(">fake_contig_4".to_string(), "GGGGGG".to_string());
        let fake_contig_5 = Contig::new_contig(">fake_contig_5".to_string(), "AAATTT".to_string());
        let fake_contig_6 = Contig::new_contig(">fake_contig_6".to_string(), "TTTCCC".to_string());
        let fake_contig_7 = Contig::new_contig(">fake_contig_7".to_string(), "CCCGGG".to_string());
        let fake_contig_8 = Contig::new_contig(">fake_contig_8".to_string(), "GGGAAA".to_string());
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


        let mut fake_bin = Bin {
        bin_contigs: Vec::new(),
        completeness: 100.0,
        contamination: 8.75,
        bin_type: BinType::prokaryote,
        bin_hash: "fake_bin_1".to_string(),

        };
        fake_bin
    }
    fn create_fake_bin2_for_unit_test() -> Bin {


        let mut fake_bin2 = Bin {
        bin_contigs: Vec::new(),
        completeness: 100.0,
        contamination: 9.0,
        bin_type: BinType::prokaryote,
        bin_hash: "fake_bin_2".to_string(),


        };
        fake_bin2
    }
    fn create_fake_bin3_for_unit_test() -> Bin {
        let mut fake_bin3 = Bin {
        bin_contigs: Vec::new(),
        completeness: 40.0,
        contamination: 9.0,
        bin_type: BinType::prokaryote,
        bin_hash: "fake_bin_3".to_string(),


        };
        fake_bin3
    }
    fn create_fake_bin4_for_unit_test() -> Bin {
        let mut fake_bin4 = Bin {
        bin_contigs: Vec::new(),
        completeness: 40.0,
        contamination: 9.0,
        bin_type: BinType::prokaryote,
        bin_hash: "fake_bin_4".to_string(),

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
        bin_1.bin_contigs = fake_contigs_for_bin_1;
        bin_2.bin_contigs = fake_contigs_for_bin_2;
        bin_3.bin_contigs = fake_contigs_for_bin_3;
        (vec![bin_1, bin_2, bin_3], fake_contigs)
    }


        #[test]
        fn test_jaccard_similarity() {
            let (fake_bins, _) = create_bin_contig_mock();
            let bin_one_against_two_result = 0.0;
            let bin_one_against_three_result = 0.125;
            let bin_two_against_three_result = 0.8;
            assert_eq!(ClusteringPrep::calc_jaccard_distance(&fake_bins[0], &fake_bins[1]), bin_one_against_two_result);
            assert_eq!(ClusteringPrep::calc_jaccard_distance(&fake_bins[0], &fake_bins[2]), bin_one_against_three_result);
            assert_eq!(ClusteringPrep::calc_jaccard_distance(&fake_bins[1], &fake_bins[2]), bin_two_against_three_result);
            let connected_bins = ClusteringPrep::get_bin_pairs_with_less_than_max_jaccard_distance(&fake_bins, 0.7);
            assert_eq!(connected_bins.len(), 1);
            assert!(connected_bins[0] == (&fake_bins[1], &fake_bins[2]));

        }
        #[test]
        fn test_graph_creation() {
            let (fake_bins, _) = create_bin_contig_mock();
            let connected_bins = ClusteringPrep::get_bin_pairs_with_less_than_max_jaccard_distance(&fake_bins, 0.7);
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
        /* 
        #[test]
        fn real_example_bins(){
            
            let bin_generator = BinGen::initialise(Some(CHECKM2_DB_PATH.to_path_buf()), Some(COMPLEASM_DB_LIB.to_path_buf()), TEST_DATA_HASH.clone(), 100.0, 0.0);
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
        */

        #[test]
        fn test_calculate_euclidean_distance_between_bin_pair() {
            let (mock_bins, mock_contigs) = create_bin_contig_mock();
            let contig_kmer_dict = ClusteringPrep::create_contig_kmer_dict_from_bins(3, mock_contigs);
   //         println!("{:?}", contig_kmer_dict);
  //      let test_rs = ClusteringPrep::calculate_euclidean_distance_between_bin_pair(&mock_bins[0], &mock_bins[1], Arc::new(contig_kmer_dict));
  //      println!("{}", test_rs);
        }

}
