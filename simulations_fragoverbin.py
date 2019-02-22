#!/usr/bin/env python3
# coding=utf-8

import numpy as np
from datetime import datetime
import os
import sys

sp = "mouse"  # sys.argv[1] "human"
genome = "mm10"  # sys.argv[2] "hg38"

# Dictionary of restrictions fragments
dic_frag = {}
with open("../data/"+sp+"/Digest_"+genome+"_HindIII_None.txt.rmap", 'r') as f1:
    for i in f1.readlines()[1:]:
        i = i.split("\t")
        frag = (int(i[1]), int(i[2]))
        chr = i[0]
        if chr in dic_frag.keys():
            dic_frag[chr].append(frag)
        else:
            dic_frag[chr] = [frag]

for i in dic_frag.keys():
    dic_frag[i] = list(set(tuple(x) for x in dic_frag[i]))
    dic_frag[i].sort(key=lambda x: x[0])

# Dictionary of interactions distances by fragment
dist = {}
vect_dist = []
trans = tot = 0
short_dist = []
super_dist = []

with open("../data/"+sp+"/all_interactions/all_interactions.txt") as f3:
    for i in f3.readlines()[1:]:
        i = i.split("\t")
        midbait = ((int(i[1]) + int(i[2])) / 2)
        midcontact = ((int(i[4]) + int(i[5])) / 2)
        frag = (str(i[0]) + ":" + str(i[1]) + ":" + str(i[2]))
        tot += 1
        if i[0] == i[3]:
            if 25000 < abs(midbait-midcontact) < 10000000:
                vect_dist.append(midbait - midcontact)
                if frag in dist.keys():
                    dist[frag].append(midbait - midcontact)
                else:
                    dist[frag] = [midbait - midcontact]

            else:
                if abs(midbait - midcontact) < 25000:
                    short_dist.append(frag + ':' + str(i[1]) + ':' + str(i[6]))
                if abs(midbait - midcontact) > 10000000:
                    super_dist.append(frag + ':' + str(i[1]) + ':' + str(i[6]))

        else:
            trans += 1

print('Dictionary of interactions distances by fragment : done ! ')
print(min(vect_dist), max(vect_dist))
print(np.median(vect_dist))
print("total :", tot)
print("ok:", len(vect_dist))
print('trans:', trans)
print('super dist:', len(super_dist))
print('short dist:', len(short_dist))
print("Number of fragments:", len(dist))

# Distribution of distance probabilities (bins)
distances = [int(i/5000) for i in vect_dist]     # 5000 = 2000*2 bin of 5kB
bin_proba = {}
for i in range(min(distances), max(distances)+1):
    bin_proba[i] = float(distances.count(i)/len(distances))

print("Starting simulations")
startTime = datetime.now()

## Simulations by fragment
output = open("../data/"+sp+"/simulations_10Mb_bin5kb_fragoverbin.txt", 'w')
if os.stat("../data/"+sp+"/simulations_10Mb_bin5kb_fragoverbin.txt").st_size == 0:
    output.write("chr_bait\tstart_bait\tend_bait\tstart\tend\n")

contact_simul = {}
running = 0
for bait in dist.keys():
    chr = bait.split(':')[0]
    fragment = str(bait.split(':')[0]) + '\t' + str(bait.split(':')[1]) + '\t' + str(bait.split(':')[2])
    midbait = ((int(bait.split(':')[1])+int(bait.split(':')[2])) / 2)
    possible_frag = []

    # Find all possible contacted fragments
    i = 0
    mid_frag = ((int(dic_frag[chr][i][0]) + int(dic_frag[chr][i][1])) / 2)

    while i < len(dic_frag[chr]) and midbait - 10000000 > ((int(dic_frag[chr][i][0]) + int(dic_frag[chr][i][1])) / 2):
        i += 1

    while i < len(dic_frag[chr]) and midbait - 10000000 < ((int(dic_frag[chr][i][0]) + int(dic_frag[chr][i][1])) / 2) < midbait - 25000:
        possible_frag.append(dic_frag[chr][i])
        i += 1

    while i < len(dic_frag[chr]) and midbait + 25000 > ((int(dic_frag[chr][i][0]) + int(dic_frag[chr][i][1])) / 2):
        i += 1

    while i < len(dic_frag[chr]) and midbait + 10000000 > ((int(dic_frag[chr][i][0]) + int(dic_frag[chr][i][1])) / 2):
        possible_frag.append(dic_frag[chr][i])
        i += 1

    # Associate probability to each possible frag
    dic_prob = {}
    for frag in possible_frag:
        first_bin = int((frag[0] - midbait)/5000)
        last_bin = int((frag[1] - midbait)/5000)
        p = []
        for i in range(first_bin, last_bin+1):
            if i in bin_proba.keys():
                p.append(bin_proba[i])
            else:
                p.append(float(0))

        dic_prob[str(frag[0])+'\t'+str(frag[1])] = np.mean(p)

    #print(dic_prob)

    # Adjust probability to sum to 1
    frag_prob = {frag: proba / sum(dic_prob.values()) for frag, proba in dic_prob.items()}

    # Selection of same number of unique restriction frag as in observed data
    select = np.random.choice(list(frag_prob), len(dist[bait]), p=list(frag_prob.values()), replace=False)
    for frag_simul in select:
        output.write(fragment + '\t' + frag_simul + '\n')

    running += 1
    print("Simulated fragments :", running, "of", len(dist))

print("Temps d'exécution :", datetime.now() - startTime)
output.close()

