#!/usr/bin/env python3
# coding=utf-8

import os

# Dictionnaire frag_position
dic_frag = {}
with open("../data/mouse/Digest_mm10_HindIII_None.txt.rmap",'r') as f2:
    for i in f2.readlines():
        i = i.split("\t")
        frag = (int(i[1]), int(i[2]))
        chr = i[0]
        if chr in dic_frag.keys():
            dic_frag[chr].append(frag)
        else:
            dic_frag[chr] = [frag]

# Trier les tuples pour chaque clé
for i in dic_frag.keys():
    dic_frag[i] = list(set(tuple(x) for x in dic_frag[i]))
    dic_frag[i].sort(key=lambda x: x[0])

# List bait
dic_bait = {}
with open("../data/mouse/Digest_mm10_HindIII_None.txt.baitmap",'r') as f1:
    for i in f1.readlines():
        i = i.split("\t")
        bait = str(i[0]) + "\t" + str(i[1]) + "\t" + str(i[2])
        if bait not in dic_bait.keys():
            dic_bait[bait] = 1

# Dictionnaire TSS
dic_TSS = {}
with open("../data/mouse/mouse_mm10_TSS.txt", 'r') as f3:
    for i in f3.readlines()[1:]:
        i = i.split("\t")
        chr = i[4]
        TSS = (int(i[3]), str(i[1]).strip("\n"))
        if chr in dic_TSS.keys():
            dic_TSS[chr].append(TSS)
        else:
            dic_TSS[chr] = [TSS]

for i in dic_TSS.keys():
    dic_TSS[i] = list(set(tuple(x) for x in dic_TSS[i]))
    dic_TSS[i].sort(key=lambda x: x[0])


# Overlap frag_position-TSS
dic_output = {}
for chr in dic_frag.keys():
    first_i = 0
    for frag in dic_frag[chr]:
        frag_start = frag[0]
        frag_end = frag[1]
        if chr in dic_TSS.keys():
            i = first_i
            frag = str(chr) + "\t" + str(frag_start) + "\t" + str(frag_end)

            # Initialization of first possible overlapping TSS
            while i < len(dic_TSS[chr]) and dic_TSS[chr][i][0] < frag_start:
                i += 1
            first_i = i

            # Adding all overlapping TSS with frag
            while i < len(dic_TSS[chr]) and dic_TSS[chr][i][0] <= frag_end:
                if frag in dic_output.keys():
                    dic_output[frag].append(dic_TSS[chr][i][1])
                else:
                    dic_output[frag] = [dic_TSS[chr][i][1]]
                i += 1

            # Adding fragment without overlapping TSS in last genome version
            if frag not in dic_output.keys():
                dic_output[frag] = ["no_TSS"]


output = open("../data/mouse/annot_frag_mm10.txt", 'w')
if os.stat("../data/mouse/annot_frag_mm10.txt").st_size == 0:
    output.write("chr\tstart\tend\tTSS\tbait\n")

for frag, TSS in dic_output.items():
    output.write(frag + "\t")
    for i in TSS:
        output.write(str(i) + ",")

    if frag in dic_bait.keys():
        output.write("\t" + "baited")
    else:
        output.write("\t" + "no_baited")

    output.write("\n")
