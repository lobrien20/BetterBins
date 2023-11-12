use std::{sync::Arc, collections::{HashMap, HashSet}};

use crate::{bin_info_storing::BinSet, contigs::Contig};


pub struct NtAnalyser;


impl NtAnalyser {
    // score nucleotide frequency similarity

    fn calculate_nucleotide_info_scores(&self, bin_set: BinSet, contigs: Vec<Arc<Contig>>) {
        let mut relevant_frequencies = vec![3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,25,30,35,40,45,50];
        for kmer_size in relevant_frequencies {
        
            let hashmap_of_contig_kmer_hashmaps = self.get_kmer_dict_for_contig(kmer_size, &contigs);
            let average_jaccard_distance = self.calculate_average_of_chosen_distance_for_binset(&bin_set, &hashmap_of_contig_kmer_hashmaps, calculate_kmer_jaccard_distance_between_contigs);
            let average_euclidean_distance = self.calculate_average_of_chosen_distance_for_binset(&bin_set, &hashmap_of_contig_kmer_hashmaps, calculate_kmer_euclidean_distance_between_contigs);
        
        }


    }
    
    fn get_kmer_dict_for_contig(&self, kmer_size: usize, contigs: &Vec<Arc<Contig>>) -> HashMap<Arc<Contig>, HashMap<String, i32>> {
        let mut hashmap_of_contig_kmer_hashmaps = HashMap::new();
        
        for contig in contigs {
            hashmap_of_contig_kmer_hashmaps.insert(contig.clone(), contig.get_kmer_frequences(kmer_size));
        
        }
        hashmap_of_contig_kmer_hashmaps

    }

    fn calculate_distances(bin_set: BinSet, hashmap_of_contig_kmer_hashmaps: HashMap<Arc<Contig>, HashMap<String, i32>>) {
        
    }

    fn calculate_average_of_chosen_distance_for_binset<F>(&self, bin_set: &BinSet, hashmap_of_contig_kmer_hashmaps: &HashMap<Arc<Contig>, HashMap<String, i32>>, calculate_distance_function: F) -> f64 
    where
        F: Fn(&HashMap<String, i32>, &HashMap<String, i32>) -> f64,

    {
        
        let mut bin_average_dists = Vec::new();
        
        for bin in &bin_set.best_bins {
            
            let vec_of_contig_kmer_hashmap: Vec<&HashMap<String, i32>> = bin.0.clone().into_iter().map(|contig| hashmap_of_contig_kmer_hashmaps.get(&contig).clone().unwrap()).collect();
            let average_dist = self.calculate_average_distance_of_bin(vec_of_contig_kmer_hashmap, &calculate_distance_function);
            bin_average_dists.push(average_dist);
        
        
        }
        let total_bins = bin_average_dists.len() as f64;
        let bin_set_average_dist = bin_average_dists.into_iter().sum::<f64>() / total_bins;
        bin_set_average_dist

        
    }
    fn calculate_average_distance_of_bin<F>(&self, vec_of_contig_kmer_hashmaps: Vec<&HashMap<String, i32>>, calculate_distance_function: F) -> f64 
    where
        F: Fn(&HashMap<String, i32>, &HashMap<String, i32>) -> f64,

    {
        let mut used_contig_hashmaps = Vec::new();
        let mut produced_distances: Vec<f64> = Vec::new();
        
        for contig_hashmap in &vec_of_contig_kmer_hashmaps {
            
            for contig_hashmap2 in &vec_of_contig_kmer_hashmaps {
                
                if used_contig_hashmaps.contains(&contig_hashmap2) {
                
                    continue
                }
            
                produced_distances.push(calculate_distance_function(contig_hashmap, contig_hashmap2));
            
            }
            used_contig_hashmaps.push(contig_hashmap);
        }
        
        let total_distances = produced_distances.len() as f64;
        let average_jaccard_distance: f64 = produced_distances.into_iter().sum::<f64>() / total_distances;
        average_jaccard_distance
    }

 
    

    // similarity in coverage
    // similarity on contig lengths
    // predicted eukaryotic/prokaryotic presence
    // overall completeness
    // overall contamination
    // proportion of contigs used
    // average bin size
    // average bin used


}

fn calculate_kmer_jaccard_distance_between_contigs(contig_1_kmer_frequencies: &HashMap<String, i32>,  contig_2_kmer_frequencies: &HashMap<String, i32> ) -> f64 {
        
    let contig_1_kmers: HashSet<String> = contig_1_kmer_frequencies.clone().into_keys().collect();
    let contig_2_kmers: HashSet<String> = contig_2_kmer_frequencies.clone().into_keys().collect();
    let intersection_size = contig_1_kmers.intersection(&contig_2_kmers).count();
    let union_size = contig_1_kmers.union(&contig_2_kmers).count();
    
    if union_size == 0 {
        return 0.0;
    }
    
    intersection_size as f64 / union_size as f64

}


fn calculate_nucleotide_distance_variance(contig_1_kmer_frequencies: &HashMap<String, i32>, contig_2_kmer_frequencies: &HashMap<String, i32>) -> f64 {
    let all_possible_kmers: Vec<String> = contig_1_kmer_frequencies.clone().into_keys().chain(contig_2_kmer_frequencies.clone().into_keys()).collect();
    let mut sum_of_contig_difference_squared = 0;
    
    for kmer in all_possible_kmers {
    
        let contig_frequency_tuple = (*contig_1_kmer_frequencies.get(&kmer).unwrap_or(&0), contig_2_kmer_frequencies.get(&kmer).unwrap_or(&0));
        let contig_difference_squared = (contig_frequency_tuple.0 - contig_frequency_tuple.1) ^ 2;
        let kmer_mean = (contig_frequency_tuple.0 + contig_frequency_tuple.1) / 2;
        let the_variance = (contig_frequency_tuple.0 - kmer_mean) ^ 2 + (contig_frequency_tuple.1 - kmer_mean) ^ 2;
        sum_of_contig_difference_squared += contig_difference_squared;
    
    }


    let variance = sum_of_contig_difference_squared / 2;
    variance as f64

}

fn calculate_nucleotide_distance_variance_at_bin_level(vec_of_contig_kmer_hashmaps: Vec<&HashMap<String, i32>>) -> f64 {
    let all_kmers: HashSet<&String> = vec_of_contig_kmer_hashmaps.iter().flat_map(|kmer_hashmap| kmer_hashmap.keys().collect::<Vec<&String>>()).collect();
    let mut sum_of_contig_difference_squared = 0;
    let number_of_contigs = vec_of_contig_kmer_hashmaps.len();
    
    
    for kmer in all_kmers {
        
        let contig_frequencies_of_kmer: Vec<i32> = vec_of_contig_kmer_hashmaps.iter()
            .map(|kmer_hashmap| *kmer_hashmap.get(kmer).unwrap_or(&0)).collect();
        
        let mean_kmer_frequency: f64 = (contig_frequencies_of_kmer.iter()
            .map(|freq| *freq as f64).sum::<f64>()) / number_of_contigs as f64;
    }
    0.0

}





fn calculate_kmer_euclidean_distance_between_contigs(contig_1_kmer_frequencies: &HashMap<String, i32>, contig_2_kmer_frequencies: &HashMap<String, i32>) -> f64 {
        
    let all_possible_kmers: Vec<String> = contig_1_kmer_frequencies.clone().into_keys().chain(contig_2_kmer_frequencies.clone().into_keys()).collect();
    let mut sum_of_contig_difference_squared = 0.0;
    let normalised_contig_1_kmer_frequencies = euclidean_normalisation(contig_1_kmer_frequencies);
    let normalised_contig_2_kmer_frequencies = euclidean_normalisation(contig_2_kmer_frequencies);
    
    for kmer in all_possible_kmers {
    
        let contig_frequency_tuple = (*normalised_contig_1_kmer_frequencies.get(&kmer).unwrap_or(&0.0), normalised_contig_1_kmer_frequencies.get(&kmer).unwrap_or(&0.0));
        let contig_difference_squared = (contig_frequency_tuple.0 - contig_frequency_tuple.1) * (contig_frequency_tuple.0 - contig_frequency_tuple.1);
      //  let kmer_variance = (contig_frequency_tuple.0 + contig_frequency_tuple.1) / 2; not sure yet..
    //    let normalised_contig_difference_squared = 
        sum_of_contig_difference_squared += contig_difference_squared;

    
    }
    
    let euclidean_distance = (sum_of_contig_difference_squared as f64).sqrt();
    euclidean_distance
}

fn euclidean_normalisation(contig_kmer_hashmap: &HashMap<String, i32>) -> HashMap<String, f64> {
    let mut normalised_kmer_hashmap = HashMap::new();

    let euclidean_norm = contig_kmer_hashmap.values()
    .map(|val| (val ^ 2) as f64).collect::<Vec<f64>>()
    .iter().sum::<f64>().sqrt();

    for (kmer, val) in contig_kmer_hashmap {
        let normalised_kmer_val = *val as f64 / euclidean_norm;
        normalised_kmer_hashmap.insert(kmer.clone(), normalised_kmer_val);

    }
    normalised_kmer_hashmap

}

