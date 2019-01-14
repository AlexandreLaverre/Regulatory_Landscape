#!/usr/bin/env python3
# coding=utf-8

# Dictionnaire bait_position
dic_bait = {}
with open("../data/human/Digest_hg38_HindIII_None.txt.baitmap",'r') as f1:
    for i in f1.readlines():
        i = i.split("\t")
        bait = (int(i[1]), int(i[2]), int(i[3]))            # tuple de 3
        chr = i[0]
        if chr in dic_bait.keys():
            dic_bait[chr].append(bait)                      # dictionnaire : ajout d'un tuple à une clé existante
        else:
            dic_bait[chr] = [bait]                          # dictionnaire : initiation d'une list de tuple pour une clé

# Trier les tuples pour chaque clé
for i in dic_bait.keys():
    dic_bait[i] = list(set(tuple(x) for x in dic_bait[i]))
    dic_bait[i].sort(key=lambda x: x[0])

# Dictionnaire TSS
dic_TSS = {}
with open("../data/human/human_hg38_TSS.txt", 'r') as f1:
    for i in f1.readlines()[1:]:
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


# Overlap bait_position-TSS
dic_output = {}
for chr in dic_bait.keys():
    first_i = 0
    for bait in dic_bait[chr]:
        bait_start = bait[0]
        bait_end = bait[1]
        if chr in dic_TSS.keys():
            i = first_i
            frag = str(chr) + "\t" + str(bait_start) + "\t" + str(bait_end) + "\t" + str(bait[2])

            # Initialization of first possible overlapping TSS
            while i < len(dic_TSS[chr]) and dic_TSS[chr][i][0] < bait_start:
                i += 1
            first_i = i

            # Adding all overlapping TSS with frag
            while i < len(dic_TSS[chr]) and dic_TSS[chr][i][0] <= bait_end:
                if frag in dic_output.keys():
                    dic_output[frag].append(dic_TSS[chr][i][1])
                else:
                    dic_output[frag] = [dic_TSS[chr][i][1]]
                i += 1

            # Adding baited fragment without overlapping TSS in last genome version
            if frag not in dic_output.keys():
                dic_output[frag] = "no_TSS"


output = open("../data/human/annot_bait_hg38_new.txt", 'w')
for frag in dic_output.keys():
    output.write(frag + "\t" + str(dic_output[frag]) + "\n")