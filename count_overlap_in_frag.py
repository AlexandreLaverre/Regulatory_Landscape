#!/usr/bin/env python3
# coding=utf-8

import numpy as np


"""
fragment = 0
count_enh = []
with open("/home/laverre/Documents/Regulatory_Landscape/data/mouse/overlap/mouse_frag_overlap_CAGE_peaks.txt") as f3:
    for i in f3.readlines()[1:]:
        i = i.strip('\n')
        i = i.split('\t')
        enh = i[3]
        if enh != 'NA':
            enh = enh.split(',')
            count_enh.append(len(enh))
        else:
            count_enh.append(0)

print("Mean enh/frag:", np.mean(count_enh))
print("Median enh/frag:", np.median(count_enh))
"""

dic_enh = {}
inter = {}
with open("/home/laverre/Documents/Regulatory_Landscape/data/mouse/all_interactions/all_interactions.txt") as f3:
    for i in f3.readlines()[1:]:
        i = i.strip('\n')
        i = i.split('\t')
        frag = str(i[0])+':'+str(i[1])+'-'+str(i[2])
        enh = str(i[3])+':'+str(i[4])+'-'+str(i[5])
        if frag not in inter.keys():
            inter[frag] = [enh]
        else:
            inter[frag].append(enh)

        if enh not in dic_enh.keys():
            dic_enh[enh] = 1

print("Nb frag:", len(inter.keys()))
print("Nb PIR:", len(dic_enh.keys()))
nb_contact = [len(inter[x]) for x in inter.keys()]

print("Mean:", np.mean(nb_contact))
print("Med:", np.median(nb_contact))


