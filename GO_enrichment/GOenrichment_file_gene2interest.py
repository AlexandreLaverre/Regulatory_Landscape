#!/usr/bin/env python3
# coding=utf-8

import os
import numpy as np

sp = "human"
data = ""  # or "_simul"
path = "/home/laverre/Documents/Regulatory_Landscape/"
path_data = path+"data/"+sp+"/"

# Genes to Baits
gene2bait = {}
with open(path_data+"overlap/"+sp+"_bait_overlap_TSS_1Kb.txt") as f1:
    for i in f1.readlines()[1:]:
        i = i.strip("\n")
        i = i.split("\t")
        bait = i[0].strip(':+')
        for gene in i[5].split(','):
            if bait not in gene2bait.keys():
                gene2bait[gene] = [bait]
            else:
                gene2bait[gene].append(bait)

# Baits to contacted fragments
if data == "_simul":
    infile = path_data + "/Simulations/simulations_" + sp + "_10Mb_bin5kb_fragoverbin_chr_merged.txt"
else:
    infile = path_data + "/all_interactions/all_interactions_chr_merged.txt_cell_names"

bait2contact = {}
with open(infile) as f2:
    for i in f2.readlines()[1:]:
        i = i.strip("\n")
        i = i.split("\t")
        bait = str(i[0]+':'+i[1]+':'+i[2])
        contact = str(i[3]+':'+i[4]+':'+i[5])
        dist = i[7]
        interest = i[9]

        if dist != 'NA' and 25000 < float(dist) < 10000000:
            if bait not in bait2contact.keys():
                bait2contact[bait] = [interest]
            else:
                bait2contact[bait].append(interest)



# Contacted fragments to interest file
# contact2interest = {}
# with open(path+"result/conservation/Sequence_conservation/contacted_sequence_composition_"+sp+data+"_merged.txt") as f3:
#     for i in f3.readlines()[1:]:
#         i = i.strip("\n")
#         i = i.split("\t")
#         contact = str(i[0]+':'+i[1]+':'+i[2])
#         interest = i[3]
#         contact2interest[contact] = int(interest)

# Output : Genes to interest file
output = open("../../result/conservation/GOrilla/"+sp+"/"+sp+"_gene_contact_sup15.txt", 'w')

nb_contact = []
for gene in gene2bait.keys():
    count_interest = []
    for bait in gene2bait[gene]:
        if bait in bait2contact.keys():
            count_interest += bait2contact[bait]

            # for contact in bait2contact[bait]:
            #     if contact in contact2interest.keys():
            #         count_interest += contact2interest[contact]

    nb_contact.append(count_interest)
    if count_interest > 15:
        output.write(gene + '\n')

output.close()

print(np.median(nb_contact))

# Exemple
# print(gene2bait["ENSMUSG00000087782"])
# print(bait2contact['chr1:13330702:13339327'])
# for contact in bait2contact['chr1:13330702:13339327']:
#    if contact in contact2interest.keys():
#        print(contact2interest[contact])

