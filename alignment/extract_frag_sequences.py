#!/usr/bin/env python3
# coding=utf-8

from Bio import SeqIO

sp = "human" #sys.argv[2]
sp2 = "mouse" #sys.argv[1]

path_data = "/home/laverre/Documents/Regulatory_Landscape/data/"
ref_genome = path_data+sp+"/genome/"+sp+".GRCh38.dna_sm.primary_assembly.fa"
ref_dict = SeqIO.to_dict(SeqIO.parse(ref_genome, "fasta"))

target_genome = path_data+sp2+"/genome/"+sp2+".GRCm38.dna_sm.primary_assembly.fa"
target_dict = SeqIO.to_dict(SeqIO.parse(target_genome, "fasta"))


with open(path_data+sp+"/genome/"+sp+"2"+sp2+"_restriction_fragments_0.1.bed", "r") as f1:
    output_sp1 = open(path_data + sp + "/genome/"+sp+"_restriction_fragments_softmask.fa", "w")
    output_sp2 = open(path_data + sp2 + "/genome/"+sp+"2"+sp2+"_restriction_fragments_softmask.fa", "w")
    for i in f1.readlines():
        i = i.split("\t")
        frag_ref = i[3].split(':')
        chr_ref = frag_ref[0].strip("chr")

        if chr_ref in ref_dict.keys():
            sequence_ref = ref_dict[chr_ref]
            cut_sequence = sequence_ref[int(frag_ref[1]):int(frag_ref[2])]
            cut_sequence.id = str(i[3] + ":+")

            SeqIO.write(cut_sequence, output_sp1, "fasta")

        chr_target = i[0].strip("chr")

        if chr_target in target_dict.keys():
            sequence_target = ref_dict[chr_ref]
            cut_sequence = sequence_target[int(i[1]):int(i[2])]
            if i[5].strip('\n') == "-":
                cut_sequence = cut_sequence.reverse_complement()

            cut_sequence.id = str(i[3] + ":+")
            cut_sequence.description = "position in ref_genome"

            SeqIO.write(cut_sequence, output_sp2, "fasta")

output_sp1.close()
output_sp2.close()

