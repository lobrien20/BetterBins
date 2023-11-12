#!/usr/bin/envs python3
import os
from Bio import SeqIO
def main():
    test_bin_dir = "/home/leah/git_repos/HeuristicMagRefiner_dev/tests/new_tests/test_bins"
    fasta_files = [f"{test_bin_dir}/{file}" for file in os.listdir(test_bin_dir)]
    file_for_unique_contigs = "/home/leah/git_repos/HeuristicMagRefiner_dev/tests/new_tests/test_bins/all_unique_contigs.fa"
    unique_contigs = get_unique_contigs_from_fasta_files(fasta_files)
    print(f"{len(unique_contigs)} unique contigs")
    generate_unique_contig_file(file_for_unique_contigs, unique_contigs)




def get_contigs_from_fasta(file):
    records = list(SeqIO.parse(file, "fasta"))
    return records

def get_unique_contigs_from_fasta_files(fasta_files):
    unique_contigs = []
    unique_contigs_by_description = []
    for file in fasta_files:
        fasta_contigs = get_contigs_from_fasta(file)
        
        for contig in fasta_contigs:
                
                if contig.description not in unique_contigs_by_description:
                    
                    unique_contigs.append(contig)
                    unique_contigs_by_description.append(contig.description)

    return unique_contigs

def generate_unique_contig_file(file_path, unique_contigs):
     
    with open(file_path, "w") as output_handle:
        SeqIO.write(unique_contigs, output_handle, "fasta")


main()