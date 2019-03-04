#!/usr/bin/env python3
# coding=utf-8

from Bio import SeqIO
import sys

sp2 = "human" #sys.argv[1]
sp = "mouse" #sys.argv[2]

path_data = "/home/laverre/Documents/Regulatory_Landscape/data/"
ref_genome = path_data+sp2+"/genome/"+sp2+".GRCh38.dna.primary_assembly.fa"
ref_dict = SeqIO.to_dict(SeqIO.parse(ref_genome, "fasta"))

output = open(path_data+sp+"/genome/"+sp+"2"+sp2+"_restriction_fragments_0.2.fa", "w")

with open(path_data+sp+"/genome/"+sp+"_frag_to_"+sp2+"_0.2.lift", 'r') as f1:
    for i in f1.readlines():
        i = i.split("\t")
        chr = i[0].strip("chr")
        if chr in ref_dict.keys():
            sequence = ref_dict[chr]
            cut_sequence = sequence[int(i[1]):int(i[2])]
            if i[5].strip('\n') == "-":
                cut_sequence = cut_sequence.reverse_complement()

            cut_sequence.id = str(i[3] + ":+")
            #cut_sequence.description = "position in ref_genome"

            SeqIO.write(cut_sequence, output, "fasta")

output.close()
