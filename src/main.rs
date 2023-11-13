#![allow(unused_imports)]
use bin_generator::{EukRepBasedPredictor, BinTypePrediction, AssumeBinType, MinimumEukMarkerGenes};
use bin_info_storage::{BinInfoStorage, BinType};
use bin_sets::BinSet;
use clap::{Parser, ValueEnum, Subcommand};
use classic_best_bin::run_classic_best_bin;
use contigs::Contig;
use initialise_bins_and_contigs::initialise_tool_through_getting_original_bins_and_contigs;
use itertools::Itertools;
use log::{debug, error, info, trace, warn};
use utils::check_and_remove_bads_in_hash_directory;
use std::{time::SystemTime, path::PathBuf, sync::Arc, fs};
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
/* 
fn main() {
 // Make sure to add functions ensuring checkm2 database and compleasm database!
    let args = Args::parse();
    if !&args.results_directory.is_dir() {
        fs::create_dir(&args.results_directory).unwrap();
    }
    setup_logger(&args.results_directory.join("output_log.txt")).expect("Failed to setup logger");

    rayon::ThreadPoolBuilder::new().num_threads(args.threads).build_global().unwrap();
    let mut hash_directory = args.results_directory.join("hash_directory/");
    match args.hash_directory {
        Some(hash_dir) => hash_directory = hash_dir,
        None => {  {if !&hash_directory.is_dir() {
            fs::create_dir(&hash_directory).unwrap();
        }}}
    }
    check_and_remove_bads_in_hash_directory(&hash_directory);

    let bin_generator = BinGenerator::initialise_bin_gen(Some(args.checkm2_db_path), Some(args.compleasm_db_dir), hash_directory.clone(), args.max_contamination, args.min_completeness);
    let bin_generator_arc = Arc::new(bin_generator);
    let (initial_bins, contigs) = gather_initial_bins_and_contig_information(&args.directory_of_bins, Arc::clone(&bin_generator_arc));
    let mut bins_to_test = Vec::new();
    match args.clustering {
        true => {
            bins_to_test = run_graph_clustering(initial_bins.clone(), Arc::clone(&bin_generator_arc), 0.5);

         //   bins_to_test = run_bin_clustering_module_to_generate_all_potential_bins(&args.results_directory, initial_bins, contigs, Arc::clone(&bin_generator_arc));
        },
        false => {
            bins_to_test = initial_bins.into_iter().map(|x| (x.bin_contigs.unwrap(), x.bin_hash, x.completeness.unwrap(), x.contamination.unwrap())).collect_vec();
        }
    }

    let bin_scorer = BinScorer {contamination_weight: args.contamination_weight, completion_weight: args.completion_weight, scoring_approach: BinSetScore::SharedCompCont, evo_completion_scale_power: args.completion_power_weight, evo_contamination_scale_power: args.contamination_power_weight};

    match args.ensemble_approach {
        EnsembleApproach::classic => {
            let best_bins = run_classic_best_bin_module(bins_to_test, Arc::clone(&bin_generator_arc), &bin_scorer);
            best_bins.create_best_bin_dir_and_info_from_best_hashes(&hash_directory, &args.results_directory.join("best_bins_directory/"), true);

        },

        EnsembleApproach::evo => {
            let contig_sets: Vec<Vec<Arc<Contig>>> = bins_to_test.into_iter().map(|x| x.0.clone()).collect();
            let evo_alg_runner = initialise_evolutionary_alg_structs(args.perc_mutation, contig_sets, args.perc_crossover, args.permutations_per_generation, args.total_generations, args.parent_picker);
            run_evolutionary_alg(Arc::clone(&bin_generator_arc), Arc::new(bin_scorer), &args.results_directory, &hash_directory, evo_alg_runner);
        }

    }

}
*/ 
fn main() {
    let args = Cli::parse();
    let mut bin_type_predictor: Box<dyn BinTypePrediction> = Box::new(EukRepBasedPredictor{});
    match args.prediction_approach {
        
        BinTypePredictionApproach::assume_eukaryote => bin_type_predictor = Box::new(AssumeBinType {assumed_bin_type: BinType::eukaryote}),
        BinTypePredictionApproach::assume_prokaryote => bin_type_predictor = Box::new(AssumeBinType {assumed_bin_type: BinType::prokaryote}),
        BinTypePredictionApproach::eukrep_majority => bin_type_predictor = Box::new(EukRepBasedPredictor{}),
        BinTypePredictionApproach::minimum_eukaryote_markers => bin_type_predictor = Box::new(MinimumEukMarkerGenes {minimum_marker_gene_count: (args.num_of_compleasm_db_markers / 2)})
    
    };

    if !&args.results_directory.is_dir() {
        fs::create_dir(&args.results_directory).unwrap();
    }
    setup_logger(&args.results_directory.join("output_log.txt")).expect("Failed to setup logger");

    rayon::ThreadPoolBuilder::new().num_threads(args.threads).build_global().unwrap();
    let mut hash_directory = args.results_directory.join("hash_directory/");
    match args.hash_directory {
        Some(hash_dir) => hash_directory = hash_dir,
        None => {  {if !&hash_directory.is_dir() {
            fs::create_dir(&hash_directory).unwrap();
        }}}
    }
    check_and_remove_bads_in_hash_directory(&hash_directory);

    let initial_bin_info_dir_path = &args.results_directory.join("initial_bin_results_information_dir/");
    fs::create_dir(&initial_bin_info_dir_path);
    let bin_info_storage = BinInfoStorage::initialise_bin_info_storer();
    
    let (bin_generator, bins) = initialise_tool_through_getting_original_bins_and_contigs(&args.contigs_file_path, initial_bin_info_dir_path, args.checkm2_db_path, args.threads, &args.path_to_bin_dir, 
        args.compleasm_db_dir, args.num_of_compleasm_db_markers, "eukaryota_odb10".to_string(), &hash_directory, args.max_contamination, 
        args.min_completeness, bin_type_predictor, bin_info_storage);

    let bin_scorer = &bin_scoring::BinScorer { contamination_weight: args.contamination_weight, completion_weight: args.completion_weight };
    let best_bins = run_classic_best_bin(bins.into_iter().map(|bin| bin.bin_contigs).collect_vec(), Box::new(bin_generator), bin_scorer);
    let best_bin_set = BinSet::make_bin_set_from_bins_vec(best_bins);
    let best_bins_directory = &args.results_directory.join("best_bins_directory/");
    best_bin_set.create_best_bin_dir_and_info_from_best_hashes(&hash_directory, best_bins_directory, true);

}

#[derive(Debug, Parser)] // requires `derive` feature
#[command(name = "BetterBins")]
#[command(about = "BetterBins bin analyser", long_about = None)]
struct Cli {
    

    #[arg(short, long)]
    contigs_file_path: PathBuf,

    #[arg(short, long, default_value = "1")]
    threads: usize,

    #[arg(long, default_value = "255")]
    num_of_compleasm_db_markers: usize,

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

    #[arg(short, long, default_value = "0.5")]
    contamination_weight: f64,
    
    #[arg(short, long, default_value = "0.5")]
    completion_weight: f64,

    #[arg(short, long, default_value = "eukrep-majority")]
    prediction_approach: BinTypePredictionApproach,

    #[arg(long)]
    compleasm_db_dir: String,

    #[arg(long)]
    checkm2_db_path: String


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
        static ref COMPLEASM_DB_LIB: PathBuf = PathBuf::from("tests/new_tests/databases_for_testing/");
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




/* 
    #[test]
    fn full_test_skip_initial_checks() {
        let (temp_initial_bin_directory, temp_results_dir, hash_dir) = prepare_initial_data_for_full_test(true);
        let bin_generator = BinGenerator::initialise_bin_gen(Some(CHECKM2_DB_PATH.to_path_buf()), Some(COMPLEASM_DB_LIB.to_path_buf()), hash_dir.clone(), 100.0, 0.0);
        let bin_gen_arc = Arc::new(bin_generator);
        setup_logger(&&temp_results_dir.join("output_log.txt")).expect("Failed to setup logger");
        rayon::ThreadPoolBuilder::new().num_threads(5).build_global().unwrap();

        let (bins, contigs) = run_initial_bin_gen_module_test(&temp_initial_bin_directory);


    }


    fn prepare_initial_data_for_full_test(skip_initial_quality_check: bool) -> (PathBuf, PathBuf, PathBuf) {
        let temp_dir_for_test = FULL_TEST_DIR.join("temp_full_test_run/");
        let temp_results_dir = temp_dir_for_test.join("results/");
        let hash_dir = temp_dir_for_test.join("hash_results/");
        let temp_initial_bin_directory = temp_dir_for_test.join("full_test_initial_bins/");
        if skip_initial_quality_check == false {
            if temp_dir_for_test.is_dir() {
            fs::remove_dir_all(&temp_dir_for_test).unwrap();
            }
            let options = CopyOptions::new();
            fs::create_dir(&temp_dir_for_test);
            fs_extra::dir::copy(FULL_TEST_DIR.join("full_test_initial_bins/"), &temp_dir_for_test, &options).unwrap();


            fs::create_dir(&temp_results_dir);
            fs::create_dir(&hash_dir);


            fs::create_dir(&temp_initial_bin_directory);

        }





        (temp_initial_bin_directory, temp_results_dir, hash_dir)
    }
    */

    fn run_initial_bin_gen_module_test(contig_file_path: &PathBuf, bin_directory_path: &PathBuf, output_directory_path: &PathBuf, hash_directory_path: &PathBuf) -> (BinGen, Vec<Bin>) {
        let bin_info_storage = BinInfoStorage::initialise_bin_info_storer();
        let (bin_gen, bins) = initialise_tool_through_getting_original_bins_and_contigs(contig_file_path, output_directory_path, CHECKM2_DB_PATH.clone().into_os_string().into_string().unwrap(), 
            6, bin_directory_path,
            COMPLEASM_DB_LIB.clone().into_os_string().into_string().unwrap(), 
            255, "eukaryota_odb10".to_string(), hash_directory_path, 
            100.0, 0.0, Box::new(EukRepBasedPredictor{}),
            bin_info_storage);


    
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
        
        let best_bin_struct = BinSet::make_bin_set_from_bins_vec(best_bins);
        best_bin_struct.create_best_bin_dir_and_info_from_best_hashes(hash_dir, classic_best_bin_test_res_dir, true);
        let files_in_best_bin_dir: Vec<DirEntry> = classic_best_bin_test_res_dir.read_dir().unwrap().filter_map(|x| x.ok()).collect();
        
        assert_eq!(files_in_best_bin_dir.len(), 4);
        let best_bin_file_path = classic_best_bin_test_res_dir.join("best_bins_information.tsv");
        let best_bin_info_results = fs::read_to_string(best_bin_file_path).expect("Unable to read best bin file");
        assert_eq!(best_bin_info_results.lines().collect::<Vec<&str>>().len(), 4);
        
    
    }
    
    /* 

    fn run_clustering_module_test(results_directory: &PathBuf, bins: Vec<Bin>, contigs: Vec<Arc<Contig>>, arc_bin_gen: Arc<BinGenerator>) -> Vec<(Vec<Arc<Contig>>, String, f64, f64)> {
        if results_directory.join("cluster_output/").is_dir() {
            fs::remove_dir_all(&results_directory.join("cluster_output/")).unwrap();
        }
        println!("{}", results_directory.to_string_lossy());
        let generated_bins = run_bin_clustering_module_to_generate_all_potential_bins(results_directory, bins, contigs, arc_bin_gen);
        let expected_cluster_file_path = &results_directory.join("cluster_output/output_cluster_file_path.tsv");
        assert!(&expected_cluster_file_path.is_file());
        let cluster_results = fs::read_to_string(expected_cluster_file_path).expect("Unable to read file");
        let expected_filtered_cluster_results = fs::read_to_string(&results_directory.join("cluster_output/output_cluster_file_eukaryotes_and_prokaryotes_split.tsv")).expect("unable to read euk prok filtered file");
        assert!(expected_cluster_file_path.is_file());
        let mut unique_ids = Vec::new();
        for res in expected_filtered_cluster_results.lines().skip(1) {
            let res_vec: Vec<&str> = res.split("\t").collect();

            println!("{}", &res);
            if !unique_ids.contains(&res_vec[1]) {
                unique_ids.push(&res_vec[1]);
            }
        }
        assert_eq!(unique_ids.len(), 3);
        generated_bins
    }


    fn run_graphing_module_test(results_directory: &PathBuf, bins: Vec<Bin>, contigs: Vec<Arc<Contig>>, arc_bin_gen: Arc<BinGenerator>) -> Vec<(Vec<Arc<Contig>>, String, f64, f64)>{
       let produced_bins = run_graph_clustering(bins, arc_bin_gen, 0.5);
       produced_bins
    }
    fn run_best_bin_classic_test(contig_sets: Vec<(Vec<Arc<Contig>>, String, f64, f64)>, bin_gen_arc: Arc<BinGenerator>, hash_directory: &PathBuf, results_directory: &PathBuf) -> BinSet {

        let bin_scorer = BinScorer { contamination_weight: 0.5, completion_weight: 0.5, scoring_approach: BinSetScore::CompContWeighted, evo_completion_scale_power: 1.0, evo_contamination_scale_power: 1.0};
        let best_bins = run_classic_best_bin_module(contig_sets, bin_gen_arc, &bin_scorer);
        let total_contigs_used = 0;
        let best_bin_result_dir_path = PathBuf::from(results_directory).join("best_bins/");
        if best_bin_result_dir_path.is_dir() {
            fs::remove_dir_all(&best_bin_result_dir_path);
        }
        best_bins.create_best_bin_dir_and_info_from_best_hashes(hash_directory, &best_bin_result_dir_path, true);
        let files_in_best_bin_dir: Vec<DirEntry> = best_bin_result_dir_path.read_dir().unwrap().filter_map(|x| x.ok()).collect();
        assert_eq!(files_in_best_bin_dir.len(), 4);
        let best_bin_file_path = best_bin_result_dir_path.join("best_bins_information.tsv");
        let best_bin_info_results = fs::read_to_string(best_bin_file_path).expect("Unable to read best bin file");
        assert_eq!(best_bin_info_results.lines().collect::<Vec<&str>>().len(), 4);
        best_bins
    }


    fn create_test_evo_alg_runner(contig_sets: Vec<(Vec<Arc<Contig>>, String, f64, f64)>) -> evo_alg_gen_tracking::EvoAlgRunner {
        let test_potential_bins: Vec<Vec<Arc<Contig>>> = contig_sets.into_iter().map(|x| x.0).collect();
        initialise_evolutionary_alg_structs(5.0, test_potential_bins, 10.0, 6, 10, "nsgaII".to_string())
    }


    fn run_evo_alg_test(contig_sets: Vec<(Vec<Arc<Contig>>, String, f64, f64)>, bin_gen_arc: Arc<BinGenerator>, hash_directory: &PathBuf, results_directory: &PathBuf) {

        let mutation_val = 1.0;
        let bin_scorer_arc = Arc::new(BinScorer { contamination_weight: 0.5, completion_weight: 0.5, scoring_approach: BinSetScore::CompContWeighted, evo_completion_scale_power: 1.0, evo_contamination_scale_power: 1.0});
        let mut evo_alg_runner = create_test_evo_alg_runner(contig_sets);
        run_evolutionary_alg(bin_gen_arc, bin_scorer_arc, results_directory, hash_directory, evo_alg_runner);

    }
*/ 
}
