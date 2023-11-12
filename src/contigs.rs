use std::collections::HashMap;


#[derive(PartialEq, Clone, Hash, Eq, Debug, PartialOrd)]
pub struct Contig {
    // contig struct containing header (the name) and its sequence
    pub header: String,
    pub sequence: String,
    pub prokaryotic_contig_info: Option<ProkaryoticContigInformation>,
    pub eukaryotic_contig_info: Option<EukaryoticContigInformation>,
    pub prok_or_euk: Option<ContigType>
}

impl Contig {
    pub fn new_contig(header_add: String, sequence_add: String) -> Contig {
        // generates contig
        Contig {
            header: header_add,
            sequence: sequence_add,
            prokaryotic_contig_info: None,
            eukaryotic_contig_info: None,
            prok_or_euk: None
        }
    }
    
    pub fn add_diamond_and_protein_information(&mut self, diamond_lines: Vec<String>, predicted_protein_lines: Vec<String>) {
    
        self.prokaryotic_contig_info = Some(ProkaryoticContigInformation { diamond_lines:  diamond_lines, predicted_protein_lines: predicted_protein_lines})
    
    }
    
    pub fn add_eukaryotic_busco_info(&mut self, busco_ids: Vec<String>) {
    
        self.eukaryotic_contig_info = Some(EukaryoticContigInformation { complete_buscos: busco_ids })
    
    }

    pub fn get_kmer_frequences(&self, length_of_kmers: usize) -> HashMap<String, i32>{
        let mut kmer_hashmap = HashMap::new();
        let mut current_kmer = String::new();
        for nt in self.sequence.chars() {
            current_kmer.push(nt);
            
            if current_kmer.len() == length_of_kmers {
                match kmer_hashmap.get(&current_kmer) {
                    
                    Some(res) => {
                    
                        let new_count = res + 1;
                        kmer_hashmap.insert(current_kmer, new_count);
                    
                    }
                    
                    None => {
                    
                        kmer_hashmap.insert(current_kmer, 1);
                    
                    }
                }
                current_kmer = String::new();  
            
            }
        }
        kmer_hashmap
    
    }


}


#[derive(PartialEq, Clone, Hash, Eq, Debug, PartialOrd)]
pub struct ProkaryoticContigInformation {
    
    pub diamond_lines: Vec<String>,
    pub predicted_protein_lines: Vec<String>

}
#[derive(PartialEq, Clone, Hash, Eq, Debug, PartialOrd)]
pub struct EukaryoticContigInformation {
    
    pub complete_buscos: Vec<String>,

}


#[derive(PartialEq, Clone, Hash, Eq, Debug, PartialOrd, Copy)]
pub enum ContigType {
    Eukaryote,
    Prokaryote
}