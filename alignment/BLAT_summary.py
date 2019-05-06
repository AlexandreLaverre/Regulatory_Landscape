#!/usr/bin/env python3
# coding=utf-8

import os
import sys
import numpy as np
from Bio import SeqIO

# sys.argv[1] from : ls 'data/"+sp1+"/genome/"+sp1+"_restriction_fragments_mask/'

sp = "human"   # sys.argv[2]
path = "/home/laverre/Documents/Regulatory_Landscape/"
path_result = path + "result/BLAT_duplication/"
path_data = path + "data/"+sp+"/genome/restriction_fragments_mask/"


seq_file = path_data + "human_restriction_fragments_mask.fa_chr1"
BLAT_file = path_result + "BLAT_output/human_restriction_fragments_mask.fa_chr1_output.psl"   # path_result + "BLAT_output/" + sys.argv[1] + "_output.psl"
output_file = path_result + "human_restriction_fragments_mask.fa_chr1_blat_summary.txt2"

# Counting N in frag
seq_dict = SeqIO.to_dict(SeqIO.parse(seq_file, "fasta"))
N_dict = {}
for frag in seq_dict.keys():
    sequence = seq_dict[frag]
    nb_N = sequence.seq.count("N")
    frag = frag.strip(':+')
    frag = frag.strip('chr')
    N_dict[frag] = nb_N


# Looking for available BLAT match
dic_match = {}
dic_part_N = {}
with open(BLAT_file) as f1:
    for i in f1.readlines()[5:]:
        i = i.strip("\n")
        i = i.split("\t")
        frag = str(i[9]).strip(':+')
        frag = frag.strip('chr')
        length_without_N = int(i[10])-N_dict[frag]
        dic_part_N[frag] = str(int(i[10])) + '\t' + str(np.around(N_dict[frag]/int(i[10]), decimals=3))
        portion = np.around(int(i[0])/length_without_N, decimals=3)
        match_ID = (str(i[13])+':'+str(i[15])+':'+str(i[16]))

        if portion > 0.2:
            if frag not in dic_match:
                dic_match[frag] = [(match_ID, portion)]

            else:
                dic_match[frag].append((match_ID, portion))

# Sorting match dictionaries by score
for k in dic_match.keys():
    dic_match[k] = list(set(tuple(x) for x in dic_match[k]))
    dic_match[k].sort(key=lambda x: x[1], reverse=True)

# Writting output
output = open(output_file, 'w')
if os.stat(output_file).st_size == 0:
    output.write("frag_ID\tlength\tpart_N\tnb_match\treal_frag\tother_score\n")

for frag in dic_match.keys():
    nb_match = len(dic_match[frag])
    all_match_ID = [match[0] for match in dic_match[frag]]
    all_match_score = [match[1] for match in dic_match[frag]]

    if frag in all_match_ID:
        real_frag = "exact_frag"
    else:
        try:
            next(i for i in all_match_score if i == 1)
            real_frag = "exact_without_N"
        except StopIteration:
            real_frag = "best_frag_score:" + str(all_match_score[0])

    if nb_match > 1:
        del all_match_score[0]
        output.write(frag + "\t" + dic_part_N[frag] + "\t" + str(nb_match) + "\t" + real_frag + "\t" + str(all_match_score) + "\n")

    else:
        output.write(frag + "\t" + dic_part_N[frag] + "\t" + str(nb_match) + "\t" + real_frag + "\t" + "None" + "\n")


output.close()

