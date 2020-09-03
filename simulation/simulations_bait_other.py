#!/usr/bin/env python3
# coding=utf-8

import numpy as np
from datetime import datetime
import os
import sys

sp = sys.argv[1] #"human"
genome = sys.argv[2] #"hg38"
cell_type = sys.argv[3]


path_data = "/beegfs/data/alaverre/Regulatory_landscape/data/"
path_result = "/beegfs/data/alaverre/Regulatory_landscape/result/simulations/"

# Dictionary of restrictions fragments
dic_frag = {}
with open(path_data+"CHICAGO_files/"+sp+"/"+genome+"/Digest_"+genome+"_HindIII_None.txt.rmap", 'r') as f1:
    for i in f1.readlines():
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

bait_list = []
with open(path_data+"CHICAGO_files/"+sp+"/"+genome+"/Digest_"+genome+"_HindIII_None.txt.baitmap", 'r') as f1:
    for i in f1.readlines():
        i = i.split("\t")
        bait = str(i[0] + ':' + i[1] + ':' + i[2])
        bait_list.append(bait)

# Dictionary of interactions distances by fragment
dist = {}
vect_dist = []
trans = tot = bait_bait = 0
short_dist = []
super_dist = []

with open(path_data+"/all_interactions/"+sp+"_samples/"+sp+"_"+cell_type+".ibed") as f3:
    for i in f3.readlines()[1:]:
        i = i.split("\t")
        midbait = ((int(i[1]) + int(i[2])) / 2)
        midcontact = ((int(i[4]) + int(i[5])) / 2)
        frag = (str(i[0]) + ":" + str(i[1]) + ":" + str(i[2]))
        tot += 1
        if i[0] == i[3]:  # cis interaction
            if str(i[3] + ':' + i[4] + ':' + i[5]) not in bait_list:  # bait-other
                if 25000 < abs(midbait-midcontact) < 10000000:        # long-range interaction
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
                bait_bait += 1
        else:
            trans += 1

print("##### Statistics #####")
print("Min and Max dist:", min(vect_dist), max(vect_dist))
print("Dist median: ", np.median(vect_dist))
print("Interactions:", tot)
print("Filtered Interactions:", len(vect_dist))
print('trans interactions:', trans)
print('bait-bait interactions:', bait_bait)
print('>10Mb interactions:', len(super_dist))
print('<25Kb interactions:', len(short_dist))
print("Fragments:", len(dist))

# Distribution of distance probabilities (bins)
distances = [int(i/5000) for i in vect_dist]     # 5000 = 2000*2 bin of 5kB
bin_proba = {}
for i in range(min(distances), max(distances)+1):
    bin_proba[i] = float(distances.count(i)/len(distances))

print("##### Starting simulations #####")
startTime = datetime.now()

## Simulations by fragment
output = open(path_result+sp+"_samples/simulations_10Mb_bin5kb_fragoverbin_"+cell_type+"_bait_other.txt", 'w')
if os.stat(path_result+sp+"_samples/simulations_10Mb_bin5kb_fragoverbin_"+cell_type+"_bait_other.txt").st_size == 0:
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
    # Adjust probability to sum to 1
    frag_prob = {key: val / sum(dic_prob.values()) for key, val in dic_prob.items()}

    # Selection of same number of unique restriction frag as in observed data
    select = np.random.choice(list(frag_prob), len(dist[bait]), p=list(frag_prob.values()), replace=False)
    for frag_simul in select:
        output.write(fragment + '\t' + frag_simul + '\n')

    running += 1

    if (running/len(dist))*10 in range(11):
        print("Avancement :", (running/len(dist))*100, "%")

print("Execution time :", datetime.now() - startTime)
output.close()
