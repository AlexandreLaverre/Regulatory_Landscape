#!/usr/bin/env python3
# coding=utf-8

from Bio import SeqIO
import sys

sp = sys.argv[1]
file = sys.argv[2] # merged_interacted_fragments_obs_simul merged_samples_simulated_fragments

if sp == "human":
    genome = "hg38"
else:
    genome = "mm10"

path = "/beegfs/data/alaverre/Regulatory_landscape/test/"
path_data = path + "data/"
path_result = path + "result/" + sp + "2other/duplication_rate/" + sp + "_restriction_fragments_mask/"
ref_genome = path_data+"genomes_rm/"+genome+"_masked.fa"
ref_dict = SeqIO.to_dict(SeqIO.parse(ref_genome, "fasta"))

with open(path_data + sp + "2other/" + sp + "_" + file + ".txt", "r") as f1:
    current_chr = "chr1"
    output = open(path_result + sp + "_" + file + "_mask.fa_chr1", "w")
    for i in f1.readlines():
        i = i.split("\t")
        chr = i[0]
        if chr != current_chr:
            output.close()
            output = open(path_result + sp + "_" + file + "_mask.fa"+"_"+chr, "w")
            current_chr = chr

        if chr in ref_dict.keys():
            sequence = ref_dict[chr]
            cut_sequence = sequence[int(i[1]):int(i[2])]

        cut_sequence.id = str(i[3] + ":+")
        SeqIO.write(cut_sequence, output, "fasta")

output.close()
