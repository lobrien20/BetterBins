# BetterBins - v0.1.0

BetterBins is a bin refinement tool for use in both eukaryotes and prokaryotes. It utilises the classical greedy based bin refinement approach first utilised by DASTool with an added graph based clustering step (experimental). For prokaryotes it relies on CheckM2, and for eukaryotes it relies on Compleasm. 

Please note that this tool is still under active development and improvements and bug fixes are ongoing.


## Installation

1. Create environment using conda/mamba/micromamba
    ''' 
    micromamba create -n BetterBins -c conda-forge -c bioconda checkm2=1.0.1 compleasm=0.2.2 eukrep=0.6.7 python=3.8 
    ''' 

2. Clone this repo, ensure Cargo and Rust are both installed.

3. 
''' 
    
    cd BetterBins
    cargo install --path.

'''
4. Install the required CheckM2 databse
'''
micromamba run -n BetterBins checkm2 database --download --path *checkm2_db_here*
'''
5. Install the required compleasm database(s) using compleasm (in this case alveolata_odb10 was used, see: https://busco.ezlab.org/busco_v4_data.html)
'''
micromamba run -n BetterBins compleasm download alveolata_odb10

'''


## Running BetterBins

#Required arguments:
'''
--path-to-bin-dir
This is the path to the directory of bins (can also be a directory of bin sub directories, eg final_bins/ which contains final_bins/metabat2_bins, final_bins/MaxBin2_bins/ etc). Accepts ".fa", ".fasta", ".fna" for bin extension.

--results-directory
Directory for results to be produced.

--compleasm-db-path
Path to the compleasm database. Can be ignored if using "assume-prokaryote" prediction approach.

--checkm2-db-path
Path to the checkm2 database.Can be ignored if using "assume-eukaryote" prediction approach.
'''
#Optional arguments:
'''
-t, --threads
Number of threads to use, default is 1

-m, --max-contamination
Maximum contamination of a bin to be seen as viable, default is 999.9. Not recommended to change.

-m, --min-completeness
Minimum completeness of a bin to be seen as viable. Not recommended to change.

-c, --contamination-weight
Weight of contamination in bin scoring method, default is 0.75.

-c, --completion-weight
Weight of completeness in bin scoring method, default is 0.25.

-p, --prediction-approach
Approach to predict type of bins, default is eukrep majority. Other possible values: assume-eukaryote, assume-prokaryote, eukrep-majority, minimum-eukaryote-markers

-m, --min-marker-prediction-minimum-marker-num
For minimum-eukaryote-markers approach how many markers to use, default is 123 (based on universal eukaryota_odb10)




'''

# Predicting bin types
Part of BetterBins' benefit is that it can be ran to identify/refine eukaryotic bins and prokaryotes, either combined or separately. It has four possible prediction approaches to run.

- assume-prokaryote: This assumes that all bins are prokaryotic and will run prokaryotic CheckM2 bin quality checking as default during refinement, skips checking for eukaryotes, good for data in which eukaryotes are unlikely.

- assume-eukaryote: This assume that all bins are eukaryotic and will run Compleasm bin quality checking as default during refinement, skipping checking for prokaryotes, good for data which you may only be interested in the eukaryotes.

- eukrep-majority: Default bin type prediction, uses eukrep contig data to predict bin, half or more of eukrep predicted contigs being eukaryotic will predict bin type as eukaryotic.

- min-eukaryote-markers: Predicts bin type based on the minimum number of eukaryotic markers. Default is 123/255 (based on eukaryota_odb10, change if using a different BUSCO marker database).

# Compleasm busco datasbase

Compleasm/BUSCO have multiple databases for different lineages for specialised markers (https://busco.ezlab.org/list_of_lineages.html). The eukaryota_odb10 is a viable choice as a universal marker set for this bin refinement tool, but choosing a specialised marker database is a better choice if you know what type of eukaryotes to expect in the metagenomic data being analysed.


# Clustering (Experimental)
Tool additionally has an iterative graph-based clustering method to generate combined bins. This is still undergoing testing so use at your own risk.

'''
-r, --run-clustering
add this to run clustering stage

-m, --max-jaccard-distance
Maximum jaccard distance between two bins (based on shared contigs) to be tested for viability

-m, --max-euclidean-distance
Maximum euclidean distance between two bins (based on kmer frequency ratios between bins) to be tested for viability

-e, --euclidean-kmer-size
Kmer size for euclidean distance calculation, default is 4 (tetranucleotide frequency)

'''
