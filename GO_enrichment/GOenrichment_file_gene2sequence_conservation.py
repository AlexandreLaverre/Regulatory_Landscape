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
    infile = path_data + "/Simulations/simulations_" + sp + "_10Mb_bin5kb_fragoverbin_chr_merged.txt"
else:
    infile = path_data + "/all_interactions/all_interactions_chr_merged.txt"

bait2contact = {}
with open(infile) as f2:
    for i in f2.readlines()[1:]:
        i = i.strip("\n")
        i = i.split("\t")
        bait = str(i[0]+':'+i[1]+':'+i[2])
        contact = str(i[3]+':'+i[4]+':'+i[5])
        if i[6] == "unbaited":
            if bait not in bait2contact.keys():
                bait2contact[bait] = [contact]
            else:
                bait2contact[bait].append(contact)

# Contacted fragments conservation
contact2interest = {}
file = "PIR_cons_all_overlap_PECAN_"+sp+"2mouse_merged.txt_onlyconserv"
with open(path+"result/conservation/Sequence_conservation/"+sp+"/contacted_sequences/"+file) as f3:
    for i in f3.readlines()[1:]:
        i = i.strip("\n")
        i = i.split("\t")
        contact = str(i[0]+':'+i[1]+':'+i[2])
        if int(i[8]) != 0:
            conservation = int(i[6])/int(i[8])  # exclude_ungapped / all_exclude
            contact2interest[contact] = conservation


# Output : Genes to interest file
output = open("../../result/conservation/GOrilla/"+sp+"/gene_seq_0.6_conserv50%.txt", 'w')
background = open("../../result/conservation/GOrilla/"+sp+"/"+sp+"_background_genes.txt", 'w')

for gene in gene2bait.keys():
    background.write(gene + '\n')
    count_contact = count_conserv = 0
    for bait in gene2bait[gene]:
        if bait in bait2contact.keys():
            for contact in bait2contact[bait]:
                if contact in contact2interest.keys():
                    count_contact += 1
                    if contact2interest[contact] > 0.6:
                        count_conserv += 1

    if count_contact != 0:
        if count_conserv/count_contact > 0.5:
            output.write(gene + '\n')

output.close()
background.close()
# Exemple
# print(gene2bait["ENSMUSG00000087782"])
# print(bait2contact['chr1:13330702:13339327'])
# for contact in bait2contact['chr1:13330702:13339327']:
#    if contact in contact2interest.keys():
#        print(contact2interest[contact])

