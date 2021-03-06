#!/usr/bin/env python3
# coding=utf-8

import os
import sys
import numpy as np
from Bio import SeqIO
import collections

# sys.argv[1] from : ls 'data/"+sp1+"/genome/"+sp1+"_restriction_fragments_mask/'

sp = sys.argv[2]
path = "/home/laverre/Documents/Regulatory_Landscape/"
path_result = path + "result/BLAT_duplication/"
path_data = path + "data/"+sp+"/genome/restriction_fragments_mask/"


seq_file = path_data + sys.argv[1]
BLAT_file = path_result + "BLAT_output/output_psl/" + sys.argv[1] + "_output.psl"
output_file = path_result + sys.argv[1] + "_blat_summary_seuil0.8.txt"

# Counting N in frag
seq_dict = SeqIO.to_dict(SeqIO.parse(seq_file, "fasta"))
N_dict = {}
for frag in seq_dict.keys():
    sequence = seq_dict[frag]
    nb_N = sequence.seq.count("N")
    frag = frag.strip(':+')
    N_dict[frag] = nb_N


# Looking for available BLAT match
dic_match = collections.defaultdict(dict)
dic_no_N = {}
with open(BLAT_file) as f1:
    for i in f1.readlines()[5:]:
        i = i.strip("\n")
        i = i.split("\t")
        frag = str(i[9]).strip(':+')

        length = int(i[10])
        length_without_N = length - N_dict[frag]
        dic_no_N[frag] = (str(length), str(length_without_N))   # str(np.around((N_dict[frag]/length)*100, decimals=2))
        match_part = int(i[0])  # np.around((int(i[0])/length_without_N)*100, decimals=2)

        match_ID = ('chr'+str(i[13])+':'+str(i[15])+':'+str(i[16]))
        match_chr = str(i[13])

        if match_part/length_without_N > 0.8:
            if frag not in dic_match.keys():
                dic_match[frag] = {}

            if match_chr not in dic_match[frag].keys():
                dic_match[frag][match_chr] = [[match_ID, match_part]]
            else:
                dic_match[frag][match_chr].append([match_ID, match_part])


### Testing overlap between matches
new_dic_match = {}
for frag in dic_match.keys():
    new_dic_match[frag] = []
    for chr_match in dic_match[frag].keys():
        if len(dic_match[frag][chr_match]) > 1:
            collapse = []
            init = dic_match[frag][chr_match][0]
            start = init[0].split(':')[1]
            end = init[0].split(':')[2]
            score = int(init[1])

            for i in range(1, len(dic_match[frag][chr_match])):
                new_frag = dic_match[frag][chr_match][i]
                new_start = new_frag[0].split(':')[1]
                new_end = new_frag[0].split(':')[2]
                new_score = int(new_frag[1])

                if new_end > end >= new_start >= start:
                    score = score + new_score
                    if new_end > end:
                        end = new_end

                else:
                    collapse.append([str(chr_match+':'+start+':'+end), score])
                    start = new_start
                    end = new_end
                    score = new_score

                if i == len(dic_match[frag][chr_match]) - 1:
                    collapse.append([str(chr_match+':'+start+':'+end), score])

            dic_match[frag][chr_match] = collapse

        for match in dic_match[frag][chr_match]:
            match[1] = float(int(match[1])/int(dic_no_N[frag][1]))*100

            new_dic_match[frag].append(tuple(match))


# Sorting match dictionaries by score
for frag in new_dic_match.keys():
    new_dic_match[frag] = list(set(tuple(x) for x in new_dic_match[frag]))
    new_dic_match[frag].sort(key=lambda x: x[1], reverse=True)

# Writting output
output = open(output_file, 'w')
#if os.stat(output_file).st_size == 0:
#    output.write("frag_ID\tlength\tno_N\tnb_match\tmatch\n")

for frag in new_dic_match.keys():
    nb_match = len(new_dic_match[frag])
    output.write(frag + "\t" + dic_no_N[frag][0] + "\t" + dic_no_N[frag][1] + "\t" + str(nb_match) + "\t")

    count = 0
    for i in new_dic_match[frag]:

        count += 1
        if count == len(new_dic_match[frag]):
            output.write(str(i) + "\n")
        else:
            output.write(str(i) + ',')

output.close()
print(sys.argv[1], "done !")
