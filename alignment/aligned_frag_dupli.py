#!/usr/bin/env python3
# coding=utf-8

import os

ref = "mouse"
sp = "human"
path = "/home/laverre/Documents/Regulatory_Landscape/"

duplication = {}
nb_frag_dupli = nb_frag_dupli_sex = 0
with open(path+"data/"+ref+"/"+ref+"_restriction_fragments_duplication_0.8.txt", 'r') as f1:
    for i in f1.readlines()[1:]:
        i = i.split('\t')
        frag = i[0]
        duplication[frag] = int(i[3])-1
        if duplication[frag] != 0:
            nb_frag_dupli += 1
            chr = frag.split(':')[0]
            if chr == "chrX" or chr == "chrY":
                nb_frag_dupli_sex += 1

print("Nb frag dupli:", nb_frag_dupli)
print("Nb frag dupli sex:", nb_frag_dupli_sex)


output = open(path+"result/alignments/frag_pairs_" + ref + "2" + sp + "_0.8.txt", 'w')
if os.stat(path+"result/alignments/frag_pairs_" + ref + "2" + sp + ".txt").st_size == 0:
    output.write("ID."+ref+"\tID."+sp+"\tnb.duplication."+ref+"\n")

missing_frag = total = 0
with open(path+"result/alignments/"+ref+"2"+sp+"/AlignmentStatistics_TBA_"+ref+"2"+sp+"_withoutnull.txt", 'r') as f1:
    for i in f1.readlines()[1:]:
        i = i.split('\t')
        frag = i[0].strip(':+')
        total += 1
        if frag in duplication.keys():
            output.write(i[0] + '\t' + i[1] + '\t' + str(duplication[frag]) + '\n')
        else:
            missing_frag += 1

output.close()
print("total frag:", total)
print("missing frag:", missing_frag)
