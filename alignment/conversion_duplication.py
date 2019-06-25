#!/usr/bin/env python3
# coding=utf-8

import os

sp = "human"
sp2 = "mouse"

path_data = "/home/laverre/Documents/Regulatory_Landscape/result/"
output_file = path_data+sp+"2"+sp2+'_converted_fragments_duplication.txt'


duplication = {}
with open(path_data+"BLAT_duplication/"+sp+"_restriction_fragments_duplication_0.8.txt", "r") as f1:
    for i in f1.readlines():
        i = i.strip("\n")
        i = i.split("\t")
        frag = (str(i[0])+':+')
        duplication[frag] = int(i[3])-1

output = open(output_file, 'w')
if os.stat(output_file).st_size == 0:
    output.write("ID."+sp+"\tID."+sp2+"\t"+sp+".duplication.score\n")

not_found = total = 0
with open(path_data+"alignments/"+sp+"2"+sp2+"/AlignmentStatistics_TBA_"+sp+"2"+sp2+"_withoutnull.txt", "r") as f1:
    for i in f1.readlines()[1:]:
        total += 1
        i = i.strip("\n")
        i = i.split("\t")
        if i[0] in duplication.keys():
            output.write(str(i[0])+'\t'+str(i[1])+'\t'+str(duplication[i[0]])+'\n')
        else:
            not_found += 1

print("Nb converted frag:", total)
print("Nb frag not found:", not_found)
output.close()


