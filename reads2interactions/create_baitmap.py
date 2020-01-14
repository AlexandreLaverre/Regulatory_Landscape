#!/usr/bin/env python3
# coding=utf-8

# Dictionnaire bait_position
dic_bait = {}
with open("../data/bait_position_Schoenfelder.txt",'r') as f1:
    for i in f1.readlines()[1:]:
        i = i.split("\t")
        bait = (int(i[1]), int(i[2]), i[4])
        chr = i[0].split("chr")
        chr = chr[1]
        if chr in dic_bait.keys():
            dic_bait[chr].append(bait)
        else:
            dic_bait[chr] = [bait]

for i in dic_bait.keys():
    dic_bait[i] = list(set(tuple(x) for x in dic_bait[i]))
    dic_bait[i].sort(key=lambda x: x[0])


# Dictionnaire rmap
dic_frag = {}
with open("../data/Digest_mm9_HindIII.rmap",'r') as f1:
    for i in f1.readlines()[1:]:
        i = i.split("\t")
        frag = (int(i[1]), int(i[2]), int(i[3]))
        chr = i[0]
        if chr in dic_frag.keys():
            dic_frag[chr].append(frag)
        else:
            dic_frag[chr] = [frag]

for i in dic_frag.keys():
    dic_frag[i] = list(set(tuple(x) for x in dic_frag[i]))
    dic_frag[i].sort(key=lambda x: x[0])


# Overlap rmap - bait_position
output = open("../data/auto_baitmap.txt", 'w')
for chr in dic_bait.keys():
    first_i = 0
    for bait in dic_bait[chr]:
        #print("bait:", bait)
        bait_start = bait[0]
        bait_end = bait[1]
        gene = bait[2]
        if chr in dic_frag.keys():
            i = first_i

            while i < len(dic_frag[chr]) and dic_frag[chr][i][1] < bait_start:
                # Initialisation du premier overlap possible
                i += 1

            first_i = i
            #print(dic_frag[chr][i][0],bait_end)
            while i < len(dic_frag[chr]) and dic_frag[chr][i][0] <= bait_end:
                output.write(str(chr) + "\t" + str(dic_frag[chr][i][0]) + "\t" + str(dic_frag[chr][i][1]) + "\t" + str(dic_frag[chr][i][2]) + "\t" + str(gene) + "\n")
                i += 1





