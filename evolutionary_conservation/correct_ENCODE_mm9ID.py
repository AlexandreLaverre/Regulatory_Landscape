#!/usr/bin/env python3
# coding=utf-8

import os
import numpy as np

path = "/home/laverre/Data/Regulatory_landscape/result/Supplementary_dataset3_annotations/mouse/ENCODE"
path_coord = path + "/ENCODE_mm10_enhancer_genomic_positions_mm10.bed"
path_align = "/home/laverre/Data/Regulatory_landscape/result/Supplementary_dataset6_regulatory_landscape_evolution/mouse/sequence_conservation/ENCODE/"

species = ["chicken", "rabbit", "rat", "cow", "dog", "elephant", "human", "macaque", "opossum"]

mm9_to_mm10 = {}
with open(path_coord) as f1:
    for i in f1.readlines():
        i = i.strip("\n")
        i = i.split("\t")
        mm9 = str(i[3])
        mm10 = str(i[0]) + ":" + str(i[1]) + ":" + str(i[2])
        mm9_to_mm10[mm9] = mm10

for sp in species:
    align = path_align + "AlignmentStatistics_Excluding_all_Exons_mouse2" + sp + "_ENCODE.txt"
    output = open(align + "_corrected_mm10", 'w')
    
    with open(align) as f1:
        output.write(f1.readline())

        for i in f1.readlines():
            i = i.strip("\n")
            i = i.split("\t")
            old_ID = str(i[0]).strip(":+")
            output.write(mm9_to_mm10[old_ID] + ":+" + "\t" + str('\t'.join(str(x) for x in i[1:])) + '\n')

    output.close()
