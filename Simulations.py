#!/usr/bin/env python3
# coding=utf-8

import numpy as np
from datetime import datetime
import os
import sys

sp = sys.argv[1] #"human"
genome = sys.argv[2] #"hg38"

# Dictionary of chromosomes length
f2 = open("../data/"+sp+"/chromosome_size_"+genome+".txt", 'r')
chromo_size = {}
for i in f2.readlines():
    i = i.split("\t")
    chromo_size[i[0]] = i[1].strip("\n")

# Dictionary of restrictions fragments
dic_frag = {}
with open("../data/"+sp+"/Digest_"+genome+"_HindIII_None.txt.rmap",'r') as f1:
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
            if 25000 < abs(midbait-midcontact) < 5000000:
                vect_dist.append(midbait - midcontact)
                if frag in dist.keys():
                    dist[frag].append(midbait - midcontact)
                else:
                    dist[frag] = [midbait - midcontact]

            else:
                if abs(midbait - midcontact) < 25000:
                    short_dist.append(frag + ':' + str(i[1]) + ':' + str(i[6]))
                if abs(midbait - midcontact) > 5000000:
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
pos = []
neg = []
for i in vect_dist:
    if i > 0:
        pos.append(int(i/5000))     # 5000 = 1000*2 bin of 5kB
    else:
        neg.append(int(i/5000))     # 10000 = 500*2 bin of 10kB

proba_pos = []
bin_pos = []
proba_neg = []
bin_neg = []
for i in range(min(pos), max(pos)):
    proba_pos.append(pos.count(i)/len(pos))
    bin_pos.append(i*5000)

for i in range(min(neg), max(neg)):
    proba_neg.append(neg.count(i) / len(neg))
    bin_neg.append(i*5000)


print("Starting simulations")
startTime = datetime.now()

#### Simulations by fragment
output = open("../data/"+sp+"/simulations_5Mb_bin5kb_uniq10_noreplace_test.txt", 'w')
if os.stat("../data/"+sp+"/simulations_5Mb_bin5kb_uniq10_noreplace_test.txt").st_size == 0:
    output.write("chr_bait\tstart_bait\tend_bait\tstart\tend\n")

contact_simul = {}
running = 0
missing_frag = 0

for frag in dist.keys():
    chr = frag.split(':')[0]
    fragment = str(frag.split(':')[0]) + '\t' + str(frag.split(':')[1]) + '\t' + str(frag.split(':')[2])
    midbait = ((int(frag.split(':')[1])+int(frag.split(':')[2])) / 2)
    proba_pos_bis = list(proba_pos)
    bin_pos_bis = list(bin_pos)
    i = len(proba_pos) - 1

    # Exclude positive distances exceeding chrom length
    while i >= 0:
        if midbait + bin_pos_bis[i] > int(chromo_size[chr]):
            del proba_pos_bis[i]
            del bin_pos_bis[i]
            print("Exclude pos +1 ! ")
            i = i - 1
        else:
            break

    # Exclude negatives distances out of null range
    proba_neg_bis = list(proba_neg)
    bin_neg_bis = list(bin_neg)
    i = 0
    while i < (len(proba_neg) - 1):
        if midbait + bin_neg_bis[i] < 0:
            del proba_neg_bis[i]
            del bin_neg_bis[i]
            print("Exclude neg +1 ! ")
        else:
            break

    # New proba distribution according to chrom boundaries
    bin_bis = bin_pos_bis + bin_neg_bis
    proba_bis = proba_pos_bis + proba_neg_bis
    proba_bis = [x / sum(proba_bis) for x in proba_bis]
    probs = np.array(proba_bis)

    # Random draw of distances according to proba distribution
    rand = np.random.choice(bin_bis, len(dist[frag])*10, p=probs)
    rand = rand.tolist()

    # Overlap with restriction frag
    simul = []
    for pos in rand:
        contact = midbait + int(pos)
        i = 0
        while i < (len(dic_frag[chr])-1) and dic_frag[chr][i][1] <= contact:
            i += 1

        simul.append(str(int(dic_frag[chr][i][0])) + '\t' + str(int(dic_frag[chr][i][1])))

    # Selection of same number of unique restriction frag as in observed data
    select = np.random.choice(list(set(simul)), len(dist[frag]), replace=False)
    for frag_simul in select:
        output.write(fragment + '\t' + frag_simul + '\n')

    running += 1
    print("Simulated fragments :", running, "of", len(dist))

print("Temps d'exécution :", datetime.now() - startTime)
output.close()
