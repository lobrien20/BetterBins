use std::{path::PathBuf, sync::{Arc, self, RwLock}, fs::{self, File}, thread, time::Duration, collections::HashSet, any::Any};

use log::{info, debug};

use crate::{prokaryotic_contig_gatherer::ProkaryoticBinQualityGetter, eukaryotic_contig_gatherer::EukaryoticBinQualityGetter, bin_info_storage::{BinInfoStorage, Bin, BinType, BinGenerationState}, utils::{generate_hash_from_contigs, create_bin_fasta, check_hash_directory_not_too_big}, contigs::{Contig, ContigType}};



 pub struct BinGen {
    pub prok_bin_quality_getter: Option<ProkaryoticBinQualityGetter>,
    pub euk_bin_quality_getter: Option<EukaryoticBinQualityGetter>,
    pub hash_directory: PathBuf,
    pub maximum_contamination: f64,
    pub minimum_completeness: f64,
    bin_info_storage: RwLock<BinInfoStorage>,
    pub bin_type_predictor: Box<dyn BinTypePrediction>,
    dry_run: bool

}

impl BinGen {
    pub fn initialise(prok_bin_quality_getter: Option<ProkaryoticBinQualityGetter>, euk_bin_quality_getter: Option<EukaryoticBinQualityGetter>, 
        hash_directory: PathBuf, maximum_contamination: f64, minimum_completeness: f64, bin_info_storage: BinInfoStorage, bin_type_predictor: Box<dyn BinTypePrediction>, dry_run: bool) -> BinGen {

        BinGen {prok_bin_quality_getter: prok_bin_quality_getter, euk_bin_quality_getter: euk_bin_quality_getter, 
            hash_directory: hash_directory, maximum_contamination: maximum_contamination, minimum_completeness: minimum_completeness, bin_info_storage: RwLock::new(bin_info_storage), 
            bin_type_predictor: bin_type_predictor, dry_run: dry_run}

    }






    fn create_relevant_bin_files_and_analyse_bin(&self, contigs: Vec<Arc<Contig>>, bin_hash_string: &str) -> Option<Bin> {
        if let Some(bin_type) = self.bin_type_predictor.predict_type_of_bin(&contigs) {
            
            let bin_hash_dir_path = &self.hash_directory.join(&bin_hash_string);
            let fasta_path = &bin_hash_dir_path.join(&format!("{}.fa", &bin_hash_string));
            
            fs::create_dir(&bin_hash_dir_path).unwrap();
            create_bin_fasta(&contigs, fasta_path);
            let mut bin_quality = None;
            
            match bin_type {
        
                BinType::eukaryote => bin_quality = self.euk_bin_quality_getter.as_ref().unwrap().analyse_bin(&contigs, bin_hash_dir_path),
                BinType::prokaryote => bin_quality = self.prok_bin_quality_getter.as_ref().unwrap().analyse_bin(&contigs, fasta_path, bin_hash_dir_path, &bin_hash_string)
            
            };

        

            match bin_quality {
               
                Some(bin_quality) => Some(self.construct_bin_object_from_info(contigs.clone(), bin_hash_string.to_string(), bin_type.clone(), (bin_quality.0, bin_quality.1))),      
                None => None

            }


        }
        else {

            None
        }

    

    }

    

    fn construct_bin_object_from_info(&self, contigs: Vec<Arc<Contig>>, bin_hash_string: String, bin_type: BinType, (completeness, contamination): (f64, f64)) -> Bin {
        
        if self.dry_run {
            fs::remove_dir_all(&self.hash_directory.join(&format!("{}", &bin_hash_string)));
        }

        let the_bin = Bin {
            bin_contigs: contigs, 
            completeness: completeness,
            contamination: contamination,
            bin_type: bin_type,
            bin_hash: bin_hash_string
        };
        debug!("Bin comp is {}, Bin cont is {}", the_bin.completeness, the_bin.contamination);

        the_bin
    }




    fn wait_for_other_thread_to_complete_bin(&self, bin_hash_string: &str) -> Option<Bin> {
        let mut time_out = 0;
        debug!("Waiting for bin to complete...");
        loop {
            time_out += 1;
            debug!("waiting for {} seconds...", time_out);
            if time_out == 10000 {
                debug!("Waiting for other thread to complete bin time out at {} seconds", time_out);
                panic!("Bin generator waiting for other thread - timeout!");
            }
            if let Some(bin) = self.bin_info_storage.read().unwrap().check_for_bin_via_hash(bin_hash_string) {
                debug!("Waiting for bin to complete finished - successful bin found!");
                return Some(bin);
            }
            if self.bin_info_storage.read().unwrap().check_if_failed_bin(bin_hash_string) {
                debug!("Waiting for bin to complete finished - failed bin found!");
                return None;
            }
            thread::sleep(Duration::from_millis(1000));
        }
    }


}

pub trait BinGenerator : Send + Sync {
    fn generate_new_bin_from_contigs(&self, contigs: Vec<Arc<Contig>>) -> Option<Bin>;
}


impl BinGenerator for BinGen {
    fn generate_new_bin_from_contigs(&self, contigs: Vec<Arc<Contig>>) -> Option<Bin> {
        debug!("TESTING BIN OF {} CONTIGS!", contigs.len());
        let bin_hash_string = generate_hash_from_contigs(&contigs);
        match self.bin_info_storage.write().unwrap().check_hypothetical_bin_status(&bin_hash_string) {
            BinGenerationState::InUse => {
        
                match self.wait_for_other_thread_to_complete_bin(&bin_hash_string) {
            
                    Some(bin) => return Some(bin),
                    None => return None,
            
                }}
        
            BinGenerationState::Failed => return None,
            BinGenerationState::Succeeded(bin) => return Some(bin),
            BinGenerationState::CreateBin => ()
        
        }
        debug!("New bin found!");

        let mut bin_result = None;
        match self.create_relevant_bin_files_and_analyse_bin(contigs, &bin_hash_string) {
            Some(bin) => bin_result = Some(bin),
            None => ()
        }
        
        if self.dry_run {

            fs::remove_dir_all(&self.hash_directory.join(&format!("{}", &bin_hash_string)));

        }
        
        self.bin_info_storage.write().unwrap().update_storage_based_on_bin_result(bin_result.clone(), bin_hash_string);
        


        bin_result

    }

}

pub trait BinTypePrediction: Send + Sync + Any {

    fn predict_type_of_bin(&self, contigs: &[Arc<Contig>]) -> Option<BinType>;
    fn required_contig_information(&self) -> Vec<String>;

}


pub struct EukRepBasedPredictor;



impl BinTypePrediction for EukRepBasedPredictor {
    
    fn predict_type_of_bin(&self, contigs: &[Arc<Contig>]) -> Option<BinType> {
        let contigs_with_types: Vec<ContigType> = contigs.iter().filter_map(|contig| contig.prok_or_euk).collect();
        let (num_of_eukaryote_contigs, num_of_prokaryotic_contigs) = contigs_with_types.iter()
            .fold((0, 0), |(euk, prok), contig_type| {
            match contig_type {
                ContigType::Eukaryote => (euk + 1, prok),
                ContigType::Prokaryote => (euk, prok + 1),
            }
        });
        
        info!("Bin:{} total contigs, {} unknown contigs, {} prokaryote contigs, and {} eukaryote contigs",  
            contigs.len(), (contigs.len() - contigs_with_types.len()), num_of_prokaryotic_contigs, num_of_eukaryote_contigs);
        
        if num_of_eukaryote_contigs > num_of_prokaryotic_contigs {
        
            info!("Eukaryote bin: {} total contigs, {} unknown contigs, {} prokaryote contigs, and {} eukaryote contigs",  
            contigs.len(), (contigs.len() - contigs_with_types.len()), num_of_prokaryotic_contigs, num_of_eukaryote_contigs);
        
            return Some(BinType::eukaryote)
        
        } else if num_of_eukaryote_contigs <= num_of_prokaryotic_contigs { // When there's an even number, assumes prokaryote since prokaryotes are more frequent in metagemomic data
                                                                            // as well as a eukaryote is more plausibly likely to have prokaryotic DNA
            info!("Prokaryote bin: {} total contigs, {} unknown contigs, {} prokaryote contigs, and {} eukaryote contigs",  
            contigs.len(), (contigs.len() - contigs_with_types.len()), num_of_prokaryotic_contigs, num_of_eukaryote_contigs);
            return Some(BinType::prokaryote)
        } else {
        
            println!("Unused Bin:{} total contigs, {} unknown contigs, {} prokaryote contigs, and {} eukaryote contigs",  
            contigs.len(), (contigs.len() - contigs_with_types.len()), num_of_prokaryotic_contigs, num_of_eukaryote_contigs);
            return None
        
        }
    }
    fn required_contig_information(&self) -> Vec<String> {
        vec!["contig_type_checking".to_string(), "prokaryote".to_string(), "eukaryote".to_string()]
    }
}


pub struct MinimumEukMarkerGenes {

    pub minimum_marker_gene_count: usize

}

impl BinTypePrediction for MinimumEukMarkerGenes {

    fn predict_type_of_bin(&self, contigs: &[Arc<Contig>]) -> Option<BinType> {
        
        let eukaryotic_unique_marker_genes: HashSet<String> = contigs.iter()    
            .filter_map(|contig| contig.eukaryotic_contig_info.clone())
            .map(|euk_info| euk_info.complete_buscos)
            .flatten()
            .collect();

        if eukaryotic_unique_marker_genes.len() >= self.minimum_marker_gene_count {
        
            info!("Eukaryotic Bin:{} total contigs, with {} marker genes", contigs.len(), eukaryotic_unique_marker_genes.len());  

            Some(BinType::eukaryote)
        
        } else {
            info!("Prokaryotic Bin:{} total contigs, with {} marker genes", contigs.len(), eukaryotic_unique_marker_genes.len());  

            Some(BinType::prokaryote)

        }

    }
    fn required_contig_information(&self) -> Vec<String> {
        vec!["prokaryote".to_string(), "eukaryote".to_string()]
    }

}

pub struct AssumeBinType {

    pub(crate) assumed_bin_type: BinType

}

impl BinTypePrediction for AssumeBinType {
    fn predict_type_of_bin(&self, contigs: &[Arc<Contig>]) -> Option<BinType> {
        Some(self.assumed_bin_type.clone())
    }
    fn required_contig_information(&self) -> Vec<String> {
        
        if self.assumed_bin_type == BinType::eukaryote {
            return vec!["eukaryote".to_string()]
        } 
        else {
            
            return vec!["prokaryote".to_string()]
        
        }
    }
}


#[cfg(test)]
mod tests {
    use crate::{contigs::{ProkaryoticContigInformation, ContigType, EukaryoticContigInformation}, initialise_bins_and_contigs::generate_all_contigs_from_contig_file, bin_generator::MinimumEukMarkerGenes};

    use super::*;
    use std::{fs, env, path::Path};
    use itertools::Itertools;
    use lazy_static::lazy_static;
    use regex::Regex;

    fn create_fake_contigs() -> Vec<Arc<Contig>> {
        let mut contig_vec = Vec::new();
        contig_vec.push(Arc::new(Contig {
            header: "na".to_string(),
            sequence: "na".to_string(),
            prokaryotic_contig_info: None,
            eukaryotic_contig_info: Some(EukaryoticContigInformation {complete_buscos: vec!["779909at2759".to_string()]}),
            prok_or_euk: Some(ContigType::Eukaryote),
            kmers: None
        }));
        contig_vec.push(Arc::new(Contig {
            header: "na".to_string(),
            sequence: "na".to_string(),
            prokaryotic_contig_info: None,
            eukaryotic_contig_info: Some(EukaryoticContigInformation {complete_buscos: vec!["290630at2759".to_string(), "290630at2759".to_string()]}), // duplicates so should only count once
            prok_or_euk: Some(ContigType::Eukaryote),
            kmers: None
        }));

        contig_vec.push(Arc::new(Contig {
            header: "na".to_string(),
            sequence: "na".to_string(),
            prokaryotic_contig_info: None,
            eukaryotic_contig_info: Some(EukaryoticContigInformation {complete_buscos: vec!["290630at2759".to_string(), "290630at2759".to_string()]}), // duplicates so should only count once
            prok_or_euk: Some(ContigType::Prokaryote),
            kmers: None
        }));



        contig_vec
    }
    #[test]
    fn test_bin_type_predictors() {
        let min_euk_marker_gene_tester = MinimumEukMarkerGenes {minimum_marker_gene_count: 1};
        let fake_bin_contigs = create_fake_contigs();
        let euk_marker_test_bin_type = min_euk_marker_gene_tester.predict_type_of_bin(&fake_bin_contigs).unwrap();
        assert_eq!(euk_marker_test_bin_type, BinType::eukaryote);
        let euk_rep_based_predictor = EukRepBasedPredictor;
        let eukrep_based_test = euk_rep_based_predictor.predict_type_of_bin(&fake_bin_contigs).unwrap();
        assert_eq!(eukrep_based_test, BinType::eukaryote);
        let assume_prok_bin_predictor = AssumeBinType {assumed_bin_type: BinType::prokaryote};
        let assume_prok_bin_test_bin_type = assume_prok_bin_predictor.predict_type_of_bin(&fake_bin_contigs).unwrap();
        assert_eq!(assume_prok_bin_test_bin_type, BinType::prokaryote);
    }
}