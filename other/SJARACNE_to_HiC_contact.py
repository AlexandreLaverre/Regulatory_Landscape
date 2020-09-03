#!/usr/bin/env python3
# coding=utf-8

import os
import numpy as np

enh = "ENCODE"

path = "/home/laverre/Data/Regulatory_landscape/"
HIC_obs = path + "/result/Supplementary_dataset4_genes_enhancers_contacts/human/" + enh + "/gene_" + enh + "_enhancers_original_interactions.txt"
HIC_simul = path + "/result/Supplementary_dataset4_genes_enhancers_contacts/human/" + enh + "/gene_" + enh + "_enhancers_simulated_interactions.txt"

enhancers = path + "data/activities_correlation/" + enh + "/ENCODE_contacted_pre_adipo_ID.txt"
promoters = path + "data/activities_correlation/" + enh + "/promoters_overlap_TSS_1kb.bed"


def conversion(file):
    conv = {}
    with open(file) as f1:
        for i in f1.readlines()[1:]:
            i = i.strip("\n")
            i = i.split("\t")

            if file == enhancers:
                ID_FOCS = i[3]
                ID_hg38 = str(i[0]) + ':' + str(i[1]) + ':' + str(i[2])
            else:
                ID_FOCS = i[0]

            conv[ID_FOCS] = ID_hg38 if file == enhancers else [gene for gene in i[4].split(",")]

    return conv


conv_enh = conversion(enhancers)
conv_prom = conversion(promoters)


def evaluate_proportion(file, contact_data):
    output = open(path + "result/SJARACNE/SJARACNE_to_" + contact_data + ".txt", 'w')
    output.write("promoter_ID\tenhancer_ID\tgene_ID\tenhancer_hg38_ID\n")

    HiC_dic = {}
    HiC_pairs = 0
    count = total = total_low_dist = total_sup_dist = 0
    common = []
    with open(file) as f1:
        for i in f1.readlines()[1:]:
            i = i.strip("\n")
            i = i.split("\t")
            gene = i[1]

            #if i[12] != "nan":  # 12 == pre_adipo
            HiC_pairs += 1
            if gene not in HiC_dic.keys():
                HiC_dic[gene] = [str(i[3])]
            else:
                HiC_dic[gene].append(str(i[3]))

    chromosomes = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
                   "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX"]

    for chrom in chromosomes:
        SJARACNE_links = path + "result/SJARACNE/" + chrom + "/consensus_network_ncol_prom-other.txt"

        with open(SJARACNE_links) as f2:
            for i in f2.readlines()[1:]:
                i = i.strip("\n")
                i = i.split("\t")
                prom_ID = i[3]
                enh_ID = i[2]
                total += 1
                prom_start = prom_ID.split("_")[1]
                prom_start = prom_start.split(":")[1]

                enh_start = enh_ID.split("_")[1]
                enh_start = enh_start.split(":")[1]

                midprom = (int(prom_start) + int(prom_ID.split("_")[2])) / 2
                midenh = (int(enh_start) + int(enh_ID.split("_")[2])) / 2
                dist = abs(midprom-midenh)
                if 25000 > dist:
                    total_low_dist += 1

                if dist > 10000000:
                    total_sup_dist += 1

                #print(prom_ID, "enh:", enh_ID, "dist", dist )

                if prom_ID in conv_prom.keys():
                    if enh_ID in conv_enh.keys():
                        for converted_prom in conv_prom[prom_ID]:
                            if converted_prom in HiC_dic.keys():

                                if str(conv_enh[enh_ID]) in HiC_dic[converted_prom]:
                                    output.write(prom_ID + "\t" + enh_ID + "\t")
                                    output.write(converted_prom + "\t" + conv_enh[enh_ID] + "\n")
                                    count += 1

                                    common.append(str(prom_ID)+"-"+str(enh_ID))

    output.close()
    print(contact_data, ":", count, "common pairs on", total, "SJARACNE links for :", HiC_pairs, "total HiC pairs, proportion :", round((count/HiC_pairs)*100, 3), "%")
    print("total_low_dist", total_low_dist, "total_sup_dist", total_sup_dist)
    return common


common_obs = evaluate_proportion(HIC_obs, "HIC_obs")
common_simul = evaluate_proportion(HIC_simul, "HIC_simul")

common_obs_simul = list(set(common_obs).intersection(common_simul))
print("Common SJARACNE in Simul & Obs :", len(common_obs_simul), "proportion :",
      round((len(common_obs_simul)/len(common_simul))*100, 3), "% of common simul")


