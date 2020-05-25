#!/usr/bin/env python3
# coding=utf-8

import os
import numpy as np

path = "/home/laverre/Data/Regulatory_landscape/data/activities_correlation/"

FOCS_links = path + "promoters_enhancers_links_FOCS_RoadMap.txt"
HIC_obs = path + "gene_enhancers_links_original_HiC_RoadMap_hg38.txt"
HIC_simul = path + "gene_enhancers_links_simulated_HiC_RoadMap_hg38.txt"

enh = path + "enhancershg19_to_enhancershg38_RoadMap.txt"
prom = path + "useful_promoters_to_genes.txt"


def conversion(file):
    conv = {}
    with open(file) as f1:
        for i in f1.readlines()[1:]:
            i = i.strip("\n")
            i = i.split("\t")

            if file == enh:
                ID = i[1].split(":")
                ID_hg38 = str(ID[0]) + ':' + str(int(ID[1])-1) + ':' + str(int(ID[2]))

            conv[i[0]] = ID_hg38 if file == enh else [gene for gene in i[2].split(",")]

    return conv


conv_enh = conversion(enh)
conv_prom = conversion(prom)


def evaluate_proportion(file, name):
    output = open(path + "FOCS_to_" + name + ".txt", 'w')
    output.write("promoter_ID\tenhancer_ID\tgene_ID\tenhancer_hg38_ID\n")

    HiC_dic = {}
    HiC_pairs = 0
    with open(file) as f1:
        for i in f1.readlines()[1:]:
            i = i.strip("\n")
            i = i.split("\t")
            gene = i[0]
            HiC_pairs += 1

            if gene not in HiC_dic.keys():
                HiC_dic[gene] = [str(i[1])]
            else:
                HiC_dic[gene].append(str(i[1]))

    count = 0
    with open(FOCS_links) as f2:
        for i in f2.readlines()[1:]:
            i = i.strip("\n")
            i = i.split("\t")
            prom_ID = i[0]
            enh_ID = i[1]

            if prom_ID in conv_prom.keys():
                if enh_ID in conv_enh.keys():

                    for converted_prom in conv_prom[prom_ID]:
                        if converted_prom in HiC_dic.keys():
                            if str(conv_enh[enh_ID]) in HiC_dic[converted_prom]:
                                output.write(prom_ID + "\t" + enh_ID + "\t")
                                output.write(converted_prom + "\t" + conv_enh[enh_ID] + "\n")
                                count += 1
    output.close()

    print(count, "common pairs on",  HiC_pairs, "total HiC pairs, proportion :", round((count/HiC_pairs)*100, 3), "%")


evaluate_proportion(HIC_obs, "HIC_obs")
evaluate_proportion(HIC_simul, "HIC_simul")
