#!/usr/bin/env python3
# coding=utf-8

from Bio import SeqIO
import sys

sp = "human" #sys.argv[2]
sp2 = "human" #sys.argv[1]

path_data = "/home/laverre/Documents/Regulatory_Landscape/data/"
ref_genome = path_data+sp2+"/genome/"+sp2+".GRCh38.dna_rm.primary_assembly.fa"
ref_dict = SeqIO.to_dict(SeqIO.parse(ref_genome, "fasta"))

#"+sp+"2"+sp2+"_restriction_fragments_0.4.bed", 'r') as f1:
with open(path_data+sp2+"/genome/Digest_hg38_HindIII_None.bed", "r") as f1:
    chr_old = 1
    output = open(path_data + sp2 + "/genome/human_restriction_fragments_mask.fa_chr1", "w")
    for i in f1.readlines():
        i = i.split("\t")
        chr = i[0].strip("chr")
        if chr != chr_old:
            output.close()
            output = open(path_data + sp2 + "/genome/human_restriction_fragments_mask.fa"+"_chr"+chr, "w")
            chr_old = chr
        # "+sp+"2"+sp2+"_restriction_fragments_0.4.fa", "w")

        if chr in ref_dict.keys():
            sequence = ref_dict[chr]
            cut_sequence = sequence[int(i[1]):int(i[2])]
            if i[5].strip('\n') == "-":
                cut_sequence = cut_sequence.reverse_complement()

            cut_sequence.id = str(i[3] + ":+")
            print(cut_sequence)
            #cut_sequence.description = "position in ref_genome"

            SeqIO.write(cut_sequence, output, "fasta")

output.close()
