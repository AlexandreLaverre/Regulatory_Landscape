#!/usr/bin/env python3
# coding=utf-8

from Bio import SeqIO
import sys

sp = "human" #sys.argv[1]
sp2 = "mouse" #sys.argv[2]

ref_genome = "../data/"+sp2+"/genome/"+sp2+".GRCm38.dna.primary_assembly.fa"
ref_dict = SeqIO.to_dict(SeqIO.parse(ref_genome, "fasta"))

output = open("../data/"+sp+"/genome/"+sp+"2"+sp2+"_restriction_fragments_0.1.fa", "w")

with open("../data/"+sp+"/genome/"+sp+"2"+sp2+"_restriction_fragments_0.1.bed", 'r') as f1:
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
