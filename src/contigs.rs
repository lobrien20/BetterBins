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

    pub fn get_kmers(&self, length_of_kmers: usize) -> Vec<String>{
        let mut kmer_vec = Vec::new();
        
        for (count, _) in self.sequence.chars().enumerate() {
            if count + length_of_kmers > self.sequence.len() {
                break
            }
            let kmer = &self.sequence[count..(count + length_of_kmers)];
            kmer_vec.push(kmer.to_string());

            }
        
        kmer_vec
    
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