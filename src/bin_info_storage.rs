use std::collections::HashSet;
use std::{collections::HashMap, sync::Arc};

use std::hash::{Hash, Hasher};

use itertools::Itertools;

use crate::{contigs::Contig, bin_sets::BinSet};

#[derive(Clone, Debug)]
pub struct BinInfoStorage {

    hash_id_to_bin_hashmap: HashMap<String, Bin>,
    failed_hash_ids: HashSet<String>,
    in_use_hashes: HashSet<String>
}

impl BinInfoStorage {
    pub fn initialise_bin_info_storer() -> BinInfoStorage {
        let mut bin_info_storer = BinInfoStorage {
            hash_id_to_bin_hashmap: HashMap::new(),
            failed_hash_ids: HashSet::new(),
            in_use_hashes: HashSet::new()
        };
        bin_info_storer
    }


    pub fn check_for_bin_via_hash(&self, bin_hash_string: &str ) -> Option<Bin> {

        match self.hash_id_to_bin_hashmap.get(bin_hash_string) {

            Some(bin) => return Some(bin.clone()),
            None => None

        }
        
    }

    pub fn check_hypothetical_bin_status(&mut self, bin_hash_string: &str) -> BinGenerationState {
        // checks if the bin has been created and was a 'good' bin, if it's currently being used, if its failed, or if it needs to be created
        if self.in_use_hashes.contains(bin_hash_string) {
            return BinGenerationState::InUse
        }
        
        if self.check_if_failed_bin(bin_hash_string) {
            return BinGenerationState::Failed
        }
        
        if let Some(bin) = self.check_for_bin_via_hash(bin_hash_string) {
            return BinGenerationState::Succeeded(bin)
        
        } else {
            self.put_hash_in_use(&bin_hash_string);
            return BinGenerationState::CreateBin
        
        }

    }


    pub fn add_bin_to_hashmap(&mut self, bin: Bin) {
    
        self.hash_id_to_bin_hashmap.insert(bin.bin_hash.clone(), bin);

        }

    pub fn check_if_failed_bin(&self, bin_hash_string: &str) -> bool {
        
        self.failed_hash_ids.contains(bin_hash_string) 

    }

    pub fn check_if_hash_in_use(&self, bin_hash_string: &str) -> bool {
        
        self.in_use_hashes.contains(bin_hash_string) 

    }

    pub fn put_hash_in_use(&mut self, bin_hash_string: &str) {
        self.in_use_hashes.insert(bin_hash_string.to_string());

    }
    pub fn take_hash_out_of_use(&mut self, bin_hash_string: &str) {
        self.in_use_hashes.remove(bin_hash_string);
    }
    pub fn update_storage_based_on_bin_result(&mut self, bin_result: Option<Bin>, bin_hash_string: String) {
        
        match bin_result {
        
            Some(bin) => self.add_bin_to_hashmap(bin),
            None => self.add_failed_bin_hash_to_hashset(bin_hash_string.clone())
        
        }
        
        self.take_hash_out_of_use(&bin_hash_string);
    }


    pub fn add_failed_bin_hash_to_hashset(&mut self, bin_hash_string: String) {
    
        self.failed_hash_ids.insert(bin_hash_string);
    
    }


}

pub enum BinGenerationState {
    InUse,
    Failed,
    Succeeded(Bin),
    CreateBin
}
#[derive(Debug, PartialEq, Clone, Hash)]
pub enum BinType {
    eukaryote,
    prokaryote
}

impl std::fmt::Display for BinType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            BinType::eukaryote => write!(f, "eukaryote"),
            BinType::prokaryote => write!(f, "prokaryote"),
        }
    }
}


#[derive(Clone, Debug)]
pub struct Bin {

    pub bin_contigs: Vec<Arc<Contig>>,
    pub completeness: f64,
    pub contamination: f64,
    pub bin_type: BinType,
    pub bin_hash: String,

}

impl Bin {



    pub fn get_bin_kmers_from_contig_kmer_hashmap<'a>(&self, contig_kmer_hashmap: &'a HashMap<Arc<Contig>, Vec<String>>) -> HashMap<String, i32> {
        let all_kmers: Vec<String> = self.bin_contigs.iter()
            .map(|contig| contig_kmer_hashmap.get(contig)
                    .expect("Critical Error: When running eukaryotic clustering contig not found in kmer hashmap!"))
            .flat_map(|kmers| kmers.iter().cloned())
            .collect();
        let mut kmer_hashmap = HashMap::new();
        for kmer in all_kmers {
            *kmer_hashmap.entry(kmer).or_insert(0) += 1;
        }
        kmer_hashmap
    }

    
}

// Enables use of bin in hashmap through using its contigs
impl PartialEq for Bin {
    fn eq(&self, other: &Self) -> bool {
        self.bin_hash == other.bin_hash
    }

    fn ne(&self, other: &Self) -> bool {
        !self.eq(other)
    }
}

impl Eq for Bin {}

impl Hash for Bin {

    fn hash<H: Hasher>(&self, state: &mut H) {

        self.bin_hash.hash(state);

    }

}



