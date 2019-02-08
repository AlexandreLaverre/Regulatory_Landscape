#!/usr/bin/env python3
# coding=utf-8

import os
import gzip

# Restriction fragment dictionary
dic_frag = {}
with open("../data/mouse/Digest_mm10_HindIII_None.txt.rmap",'r') as f1:
    for i in f1.readlines():
        i = i.split("\t")
        frag = (int(i[1]), int(i[2]))
        chr = i[0]
        if chr in dic_frag.keys():
            dic_frag[chr].append(frag)
        else:
            dic_frag[chr] = [frag]

# Sorting fragments
for i in dic_frag.keys():
    dic_frag[i] = list(set(tuple(x) for x in dic_frag[i]))
    dic_frag[i].sort(key=lambda x: x[0])

# Enhancer CAGE dictionary
dic_enh = {}
with gzip.open('../data/mouse/enhancers/F5.mm10.enhancers.bed.gz', 'r') as f2:
    for i in f2.readlines():
        i = i.decode('UTF-8')
        i = i.strip('\n')
        i = i.split("\t")
        chr = i[0].strip("chr")
        if chr in dic_enh.keys():
            dic_enh[chr].append([int(i[1]), int(i[2])])
        else:
            dic_enh[chr] = [[int(i[1]), int(i[2])]]

# Sorting enhancer
for i in dic_enh.keys():
    dic_enh[i].sort(key=lambda x: x[0])


# Overlap frag_position-TSS
dic_output = {}
for chr in dic_frag.keys():
    first_i = 0
    for frag in dic_frag[chr]:
        frag_start = frag[0]
        frag_end = frag[1]
        if chr in dic_enh.keys():
            i = first_i
            frag = str(chr) + "\t" + str(frag_start) + "\t" + str(frag_end)

            # Initialization of first possible overlapping enh
            while i < len(dic_enh[chr]) and dic_enh[chr][i][1] < frag_start:
                i += 1
            first_i = i

            # Adding all overlapping enh with frag
            while i < len(dic_enh[chr]) and dic_enh[chr][i][0] <= frag_end:
                M = max(frag_start, dic_enh[chr][i][0])
                m = min(frag_end, dic_enh[chr][i][1])
                if M <= m:
                    if frag in dic_output.keys():
                        dic_output[frag].append(str(dic_enh[chr][i][0])+":"+str(dic_enh[chr][i][1]))
                    else:
                        dic_output[frag] = [str(dic_enh[chr][i][0])+":"+str(dic_enh[chr][i][1])]
                i += 1

output = open("../data/mouse/overlap_frag_CAGE_enh.txt", 'w')
if os.stat("../data/mouse/overlap_frag_CAGE_enh.txt").st_size == 0:
    output.write("chr\tstart\tend\tCAGE_enh\n")

for frag, TSS in dic_output.items():
    output.write(frag + "\t")
    for i in TSS:
        output.write(str(i) + ",")
    output.write("\n")
