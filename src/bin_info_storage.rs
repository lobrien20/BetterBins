use std::collections::HashSet;
use std::{collections::HashMap, sync::Arc};

use std::hash::{Hash, Hasher};

use itertools::Itertools;

use crate::{contigs::Contig, bin_sets::BinSet};

#[derive(Clone, Debug)]
pub struct BinInfoStorage {

    hash_id_to_bin_hashmap: HashMap<String, Bin>,
    failed_hash_ids: HashSet<String>
}

impl BinInfoStorage {
    pub fn initialise_bin_info_storer() -> BinInfoStorage {
        let mut bin_info_storer = BinInfoStorage {
            hash_id_to_bin_hashmap: HashMap::new(),
            failed_hash_ids: HashSet::new()
        };
        bin_info_storer
    }


    pub fn check_for_bin_via_hash(&self, bin_hash_string: &str ) -> Option<Bin> {

        match self.hash_id_to_bin_hashmap.get(bin_hash_string) {

            Some(bin) => return Some(bin.clone()),
            None => None

        }
        
    }


    pub fn add_bin_to_hashmap(&mut self, bin: Bin) {
    
        self.hash_id_to_bin_hashmap.insert(bin.bin_hash.clone(), bin);

        }

    pub fn check_if_failed_bin(&self, bin_hash_string: &str) -> bool {
        
        self.failed_hash_ids.contains(bin_hash_string) 

    }


    pub fn add_failed_bin_hash_to_hashset(&mut self, bin_hash_string: String) {
    
        self.failed_hash_ids.insert(bin_hash_string);
    
    }


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
    pub bin_hash: String


}

impl Bin {
    pub fn get_all_bin_kmers(&self, kmer_size: usize) -> Vec<&Vec<String>> {
        self.bin_contigs.iter()
        .map(|contig| contig.kmers.as_ref().unwrap()) // unwrap the Option
        .collect()
    }
    pub fn get_bin_kmers_from_contig_kmer_hashmap(&self, contig_kmer_hashmap: &Arc<HashMap<Arc<Contig>, Vec<String>>>) -> Vec<String> {
        self.bin_contigs.iter()
        .map(|contig| contig_kmer_hashmap.get(contig).unwrap().clone())
        .flatten() // unwrap the Option
        .collect()
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



