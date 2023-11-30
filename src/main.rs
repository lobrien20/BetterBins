#![allow(unused_imports)]
use bin_generator::{EukRepBasedPredictor, BinTypePrediction, AssumeBinType, MinimumEukMarkerGenes, BinGen};
use bin_info_storage::{BinInfoStorage, BinType};
use bin_sets::BinSet;
use clap::{Parser, ValueEnum, Subcommand};
use classic_best_bin::run_classic_best_bin;
use contigs::Contig;
use graphing::run_graph_clustering;
use initialise_bins_and_contigs::initialise_tool_through_getting_original_bins_and_contigs;
use itertools::Itertools;
use log::{debug, error, info, trace, warn};
use utils::check_and_remove_bads_in_hash_directory;
use std::{time::SystemTime, path::PathBuf, sync::Arc, fs, hash};

use crate::graphing::run_additional_eukaryotic_clustering_stage;
pub mod utils;
pub mod contigs;
pub mod classic_best_bin;
pub mod eukaryotic_contig_gatherer;
pub mod initialise_bins_and_contigs;
pub mod prokaryotic_contig_gatherer;
pub mod contig_type_predictor;
pub mod bin_info_storage;
pub mod bin_sets;
pub mod bin_generator;
pub mod bin_scoring;
pub mod graphing;

fn main() {
    let args = Cli::parse();
    let mut bin_type_predictor: Box<dyn BinTypePrediction> = Box::new(EukRepBasedPredictor{});
    match args.prediction_approach {
        
        BinTypePredictionApproach::assume_eukaryote => bin_type_predictor = Box::new(AssumeBinType {assumed_bin_type: BinType::eukaryote}),
        BinTypePredictionApproach::assume_prokaryote => bin_type_predictor = Box::new(AssumeBinType {assumed_bin_type: BinType::prokaryote}),
        BinTypePredictionApproach::eukrep_majority => bin_type_predictor = Box::new(EukRepBasedPredictor{}),
        BinTypePredictionApproach::minimum_eukaryote_markers => bin_type_predictor = Box::new(MinimumEukMarkerGenes {minimum_marker_gene_count: args.min_marker_prediction_minimum_marker_num})
    
    };

    if !&args.results_directory.is_dir() {
        fs::create_dir(&args.results_directory).unwrap();
    }
    setup_logger(&args.results_directory.join("betterbins_output.txt")).expect("Failed to setup logger");

    rayon::ThreadPoolBuilder::new().num_threads(args.threads).build_global().unwrap();
    let mut hash_directory = args.results_directory.join("hash_directory/");
    match args.hash_directory {
        Some(hash_dir) => hash_directory = hash_dir,
        None => {  {if !&hash_directory.is_dir() {
            fs::create_dir(&hash_directory).unwrap();
        }}}
    }
    info!("Checking for any bad directories in hash that may have been produced in previous runs...");
    check_and_remove_bads_in_hash_directory(&hash_directory);


    let initial_bin_info_dir_path = &args.results_directory.join("initial_bin_results_information_dir/");
    fs::create_dir(&initial_bin_info_dir_path);
    let bin_info_storage = BinInfoStorage::initialise_bin_info_storer();
    let (bin_generator, bins) = initialise_tool_through_getting_original_bins_and_contigs(initial_bin_info_dir_path, args.checkm2_db_path, args.threads, &args.path_to_bin_dir, 
        args.compleasm_db_path, &hash_directory, args.max_contamination, 
        args.min_completeness, bin_type_predictor, bin_info_storage, args.dry_run);

    let bin_scorer = &bin_scoring::BinScorer { contamination_weight: args.contamination_weight, completion_weight: args.completion_weight };
    let arc_bin_gen = Arc::new(bin_generator);
    let mut bin_arc_contigs = None;
    if args.run_clustering {
       let cluster_output_directory = args.results_directory.join("cluster_results_directory");
       let bin_set = run_graph_clustering(bins, Arc::clone(&arc_bin_gen), args.max_jaccard_distance, cluster_output_directory);
       let eukaryotic_bins = bin_set.bins.iter()
        .filter(|bin| bin.bin_type == BinType::eukaryote)
        .map(|bin| Arc::clone(bin)) // Clones the Arc<Bin>
        .map(|arc_bin| (*arc_bin).clone()) // Clones the inner Bin
        .collect_vec();
        info!("Jaccard distance based clustering complete!");
       if eukaryotic_bins.len() > 0 {
			let eukaryotic_cluster_output_directory = args.results_directory.join("extra_eukaryotic_clustering_results_directory");
			let new_bin_set = run_additional_eukaryotic_clustering_stage(&eukaryotic_bins, Arc::clone(&arc_bin_gen), args.max_euclidean_distance, eukaryotic_cluster_output_directory, args.euclidean_kmer_size);
			bin_arc_contigs = Some(new_bin_set.bins.into_iter().map(|bin| bin.bin_contigs.clone()).collect_vec());

       } else {
    		bin_arc_contigs = Some(bin_set.bins.into_iter().map(|bin| bin.bin_contigs.clone()).collect_vec());
       }
       

    } else {
        bin_arc_contigs = Some(bins.into_iter().map(|bin| bin.bin_contigs).collect_vec());
    }

    let bin_generator = Arc::try_unwrap(arc_bin_gen).ok().unwrap();
    let best_bins = run_classic_best_bin(bin_arc_contigs.unwrap(), Box::new(bin_generator), bin_scorer).into_iter().map(|bin| Arc::new(bin)).collect_vec();
    let best_bin_set = BinSet::make_bin_set_from_bins_vec(best_bins);
    let best_bins_directory = &args.results_directory.join("best_bins_directory");
    if args.dry_run {
        best_bin_set.create_bin_set_dir_and_info_from_best_hashes(&hash_directory, best_bins_directory, false);

    } else {
    
        best_bin_set.create_bin_set_dir_and_info_from_best_hashes(&hash_directory, best_bins_directory, true);

    }
    fs::remove_dir_all(&hash_directory).unwrap();

}

#[derive(Debug, Parser)] // requires `derive` feature
#[command(name = "BetterBins")]
#[command(about = "BetterBins bin refinement tool", long_about = None)]
struct Cli {
    

    #[arg(short, long, default_value = "1")]
    threads: usize,

    #[arg(short, long)]
    path_to_bin_dir: PathBuf,

    #[arg(short, long)]
    results_directory: PathBuf,

    #[arg(long)]
    hash_directory: Option<PathBuf>,

    #[arg(short, long, default_value = "999.9")]
    max_contamination: f64,             
    
    #[arg(short, long, default_value = "0.0")]
    min_completeness: f64,
    
    #[arg(short, long, default_value = "false")]
    run_clustering: bool,

    #[arg(short, long, default_value = "0.75")]
    contamination_weight: f64,
    
    #[arg(short, long, default_value = "0.25")]
    completion_weight: f64,

    #[arg(short, long, default_value = "eukrep-majority")]
    prediction_approach: BinTypePredictionApproach,

    #[arg(long)]
    compleasm_db_path: String,

    #[arg(long)]
    checkm2_db_path: String,

    #[arg(short, long, default_value = "0.5")]
    max_jaccard_distance: f64,

    #[arg(short, long, default_value = "0.5")]
    max_euclidean_distance: f64,

    #[arg(short, long, default_value = "4")]
    euclidean_kmer_size: usize,

    #[arg(short, long, default_value = "123")]
    min_marker_prediction_minimum_marker_num: usize,
    
    #[arg(short, long, default_value = "false")]
    dry_run: bool




}




#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, ValueEnum, Debug)]
enum BinTypePredictionApproach {
    assume_eukaryote,
    assume_prokaryote,
    eukrep_majority,
    minimum_eukaryote_markers
}

#[derive(Debug, Subcommand)]
enum Commands {
    #[command(arg_required_else_help = true)]
    classic {

    },
    evo {
        
    }
}






#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, ValueEnum)]
enum EnsembleApproach {
    classic,
    evo
}



fn setup_logger(log_file_path: &PathBuf) -> Result<(), fern::InitError> {
    fern::Dispatch::new()
        .format(|out, message, record| {
            out.finish(format_args!(
                "[{} {} {}] {}",
                chrono::Local::now().format("[%Y-%m-%d][%H:%M:%S]"),
                record.level(),
                record.target(),
                message
            ))
        })
        .level(log::LevelFilter::Debug)
        .chain(std::io::stdout())
        .chain(fern::log_file(log_file_path)?)
        .apply()?;
    Ok(())
}

mod tests {

    use crate::{bin_info_storage::Bin, bin_generator::BinGen, bin_scoring::BinScorer};

    use super::*;
    use std::{fs::{self, DirEntry, remove_dir_all}, env::{self, temp_dir}, path::Path, hash};
    use fs_extra::dir::CopyOptions;
    use lazy_static::lazy_static;
    lazy_static! {
        static ref FULL_TEST_DIR: PathBuf = PathBuf::from("tests/new_tests/full_scale_test/");
    }
    lazy_static! { // unit test will only work if database is added here, not uploaded to git repo due to size
        static ref COMPLEASM_DB_LIB: String ="tests/new_tests/databases_for_testing/eukaryota_odb10/".to_string();
    }
    lazy_static! { // unit test will only work if database is added here, not uploaded to git repo due to size

        static ref CHECKM2_DB_PATH: PathBuf = PathBuf::from("tests/new_tests/databases_for_testing/uniref100.KO.1.dmnd");

    }
                       

    #[test]
    fn full_test() {
        // run if no initial hash directory data generated, takes longer since quality checks need to be done for all the samples. But more thoorough as tests the beginning module.

        rayon::ThreadPoolBuilder::new().num_threads(6).build_global().unwrap();

        let contig_file_path = &FULL_TEST_DIR.join("test_1/all_unique_contigs.fa");
        let bin_dir_path = &FULL_TEST_DIR.join("test_1/test_bins");
        let output_directory_path = &FULL_TEST_DIR.join("test_1/results/");
        let hash_directory_path = &FULL_TEST_DIR.join("test_1/results/hash_dir/");
        fs::create_dir(hash_directory_path);
        let (bin_gen, bins) = run_initial_bin_gen_module_test(contig_file_path, bin_dir_path, output_directory_path, hash_directory_path);
        run_classic_best_bin_module_test(bin_gen, bins, &output_directory_path.join("classic_best_bin_results/"));
    }


    fn run_initial_bin_gen_module_test(contig_file_path: &PathBuf, bin_directory_path: &PathBuf, output_directory_path: &PathBuf, hash_directory_path: &PathBuf) -> (BinGen, Vec<Bin>) {
        let bin_info_storage = BinInfoStorage::initialise_bin_info_storer();
        let (bin_gen, bins) = initialise_tool_through_getting_original_bins_and_contigs(output_directory_path, CHECKM2_DB_PATH.clone().into_os_string().into_string().unwrap(), 
            6, bin_directory_path,
            COMPLEASM_DB_LIB.to_string(), 
            hash_directory_path, 
            100.0, 0.0, Box::new(EukRepBasedPredictor{}),
            bin_info_storage, false);


    
        assert_eq!(bins.len(), 6);
        assert_eq!(bins.iter().fold(0, | acc, bin| acc + bin.bin_contigs.len()), 3181);
        (bin_gen, bins)

    }

    fn run_classic_best_bin_module_test(bin_generator: BinGen, bins: Vec<Bin>, classic_best_bin_test_res_dir: &PathBuf) {
        let contig_sets = bins.iter().map(|bin| bin.bin_contigs.clone()).collect_vec();
        let bin_scorer = BinScorer::initialise_bin_scorer(0.5, 0.5);
        let hash_dir = &bin_generator.hash_directory.clone();
        let best_bins = run_classic_best_bin(contig_sets, Box::new(bin_generator), &bin_scorer);
        println!("{}", best_bins.len());
        for bin in &best_bins { // we expect three bins of descending quality
            println!("Bin has: {} contigs, and {} completeness and {} contamination", bin.bin_contigs.len(), bin.completeness, bin.contamination);
        }
        let all_contigs = best_bins.iter().map(|bin| &bin.bin_contigs).collect_vec();
        let unique_contigs = all_contigs.iter().unique().collect_vec();
        if all_contigs.len() != unique_contigs.len() {
            panic!("classic_best_bin_module_test_error: Some bins have the same contigs!");
        }
        
        let best_bin_struct = BinSet::make_bin_set_from_bins_vec(best_bins.into_iter().map(|bin| Arc::new(bin)).collect_vec());
        best_bin_struct.create_bin_set_dir_and_info_from_best_hashes(hash_dir, classic_best_bin_test_res_dir, true);
        let files_in_best_bin_dir: Vec<DirEntry> = classic_best_bin_test_res_dir.read_dir().unwrap().filter_map(|x| x.ok()).collect();
        
        assert_eq!(files_in_best_bin_dir.len(), 4);
        let best_bin_file_path = classic_best_bin_test_res_dir.join("bin_set_information.tsv");
        let best_bin_info_results = fs::read_to_string(best_bin_file_path).expect("Unable to read best bin file");
        assert_eq!(best_bin_info_results.lines().collect::<Vec<&str>>().len(), 4);
        
    
    }
    
}
