#!/usr/bin/env python3
# coding=utf-8

sp = "human"
data = ""  # or "_simul"
file = "/home/laverre/Documents/Regulatory_Landscape/data/"+sp+"/all_interactions/all_interactions_chr_merged.txt"
output = open("/home/laverre/Documents/Regulatory_Landscape/data/"+sp+"/human_hg38all_merged_frag.txt", 'w')

with open(file, 'r') as f:
    for i in f.readlines()[1:]:
        i = i.strip("\n")
        i = i.split("\t")
        PIR = (str(i[3]) + "\t" + str(i[4]) + "\t" + str(i[5]))
        ID = i[7]

        if len(ID.split(',')) > 1:
            output.write(PIR + "\t" + ID + '\n')

output.close()
