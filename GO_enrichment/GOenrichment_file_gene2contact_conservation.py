#!/usr/bin/env python3
# coding=utf-8

import os

sp = "human"
data = ""  # or "_simul"
path = "/home/laverre/Documents/Regulatory_Landscape/"
path_data = path+"data/"+sp+"/"

# Baits to genes
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
    infile = path + "result/conservation/interaction_conservation/human_simul_to_mouse_0.1_merged_all_infos.txt_target_correction"
else:
    infile = path + "result/conservation/interaction_conservation/human_to_mouse_0.1_merged_all_infos.txt_target_correction"

bait2contact = {}
bait2contact_conserv = {}
with open(infile) as f2:
    for i in f2.readlines()[1:]:
        i = i.strip("\n")
        i = i.split("\t")
        contact = i[0]
        bait = contact.split('-')[0]

        if bait not in bait2contact.keys():
            bait2contact[bait] = 1
            bait2contact_conserv[bait] = 0
            if i[10] != "[]":
                bait2contact_conserv[bait] += 1
        else:
            bait2contact[bait] += 1
            if i[10] != "[]":
                bait2contact_conserv[bait] += 1

# Output : Genes to interest file
output = open("../../result/conservation/GOrilla/"+sp+"/gene_contact_conserv_sup25%.txt", 'w')
background = open("../../result/conservation/GOrilla/"+sp+"/"+sp+"_background_genes_contact.txt", 'w')
distrib = open("../../result/conservation/GOrilla/"+sp+"/gene_contact_distribution.txt", 'w')

for gene in gene2bait.keys():
    count_contact = count_conserv = 0
    for bait in gene2bait[gene]:
        if bait in bait2contact.keys():
            count_contact += bait2contact[bait]
            count_conserv += bait2contact_conserv[bait]

    if count_contact != 0:
        background.write(gene + '\n')
        distrib.write(gene + '\t' + str(count_conserv / count_contact) + '\n')
        if count_conserv/count_contact > 0.25:
            output.write(gene + '\n')

output.close()
background.close()

# Exemple
# print(gene2bait["ENSMUSG00000087782"])
# print(bait2contact['chr1:13330702:13339327'])
# for contact in bait2contact['chr1:13330702:13339327']:
#    if contact in contact2interest.keys():
#        print(contact2interest[contact])

