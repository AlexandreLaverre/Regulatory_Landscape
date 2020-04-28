#!/usr/bin/env python3.7
# coding=utf-8

from Bio import SeqIO
import sys
import gzip

sp = sys.argv[1] # ex : human
sp2 = sys.argv[2] # ex : opossum
file = sys.argv[3] # ex : CAGE, merged_interacted_fragments,...

path = "/beegfs/data/alaverre/Regulatory_landscape/test/data/"
path_data = path + sp + "2other/" + sp + "2" + sp2 + "/"
path_result = "/beegfs/data/alaverre/Regulatory_landscape/test/result/"+sp+"2other/"+sp+"2"+sp2+"/"

target_genome = path+"genomes_sm/"+sp2+".dna_sm.toplevel.fa.gz"
target_dict = SeqIO.to_dict(SeqIO.parse(gzip.open(target_genome, "rt"), "fasta"))

with open(path_data + sp + "2" + sp2 + "_" + file +".bed", "r") as f1:  #
    output_sp = open(path_result + sp + "2" + sp2 + "_" + file +"_softmask.fa", "w")
    for i in f1.readlines():
        i = i.split("\t")
        chr_target = i[0]

        if chr_target in target_dict.keys():
            sequence_target = target_dict[chr_target]
            cut_sequence = sequence_target[int(i[1]):int(i[2])]
            if i[5].strip('\n') == "-":
                cut_sequence = cut_sequence.reverse_complement()

            cut_sequence.id = str(i[3] + ":+")
            cut_sequence.description = "position in ref_genome"

            SeqIO.write(cut_sequence, output_sp, "fasta")

output_sp.close()
