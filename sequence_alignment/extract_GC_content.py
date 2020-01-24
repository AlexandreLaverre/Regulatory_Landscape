#!/usr/bin/env python3
# coding=utf-8

import os
import sys
from Bio import SeqIO

sp = sys.argv[1]
file = sys.argv[2]  # samples_simulations merged_interacted_fragments GRO_seq ... bait

path = "/beegfs/data/alaverre/Regulatory_landscape/test/result/" + sp + "2other/"
seq_file = path + sp + "2" + sp + "/" + sp + "2" + sp + "_" + file + "_softmask.fa"
output_file = path + file + "_stats.txt"

# Counting N in frag
seq_dict = SeqIO.to_dict(SeqIO.parse(seq_file, "fasta"))
content_dict = {}
for frag in seq_dict.keys():
    sequence = seq_dict[frag]
    nb_N = sequence.seq.count("N")
    nb_GC = sequence.seq.count("G") + sequence.seq.count("g") + sequence.seq.count("C") + sequence.seq.count("c")
    length = len(sequence)
    frag = frag.strip(':+')
    content_dict[frag] = (length, nb_N, nb_GC)

print("Writting output...")
# Writting output
output = open(output_file, 'w')
if os.stat(output_file).st_size == 0:
    output.write("frag_ID\tlength\tnb_N\tnb_GC\n")

for frag in content_dict.keys():
    output.write(
        frag + "\t" + str(content_dict[frag][0]) + "\t" + str(content_dict[frag][1]) + "\t" + str(content_dict[frag][2]) + "\n")

output.close()
