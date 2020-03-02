#!/usr/bin/env python3
# coding=utf-8

import os
import numpy as np

origin_sp = "human"
enhancers = ["GRO_seq"]
target_sps = ["chicken", "opossum", "rat", "mouse", "rabbit", "dog", "elephant", "cow", "macaque"]
datas = ["", "_simul"]

path = "/home/laverre/Documents/Regulatory_Landscape/"
path_data = path + "data/"+origin_sp+"/"
path_result = path + "result/conservation/"


def conserv_enh(origin_sp, target_sp, enhancers, data):
    print(origin_sp, "to", target_sp, "; enh data : ", enhancers, data)

    # Alignment score of each enhancer
    file = enhancers + "/AlignmentStatistics_Excluding_all_Exons_" + origin_sp + "2" + target_sp + "_" + enhancers + ".txt"
    align = {}
    with open(path_result + "alignments/" + origin_sp + "/" + file) as f2:
        for i in f2.readlines()[1:]:
            i = i.strip("\n")
            i = i.split("\t")
            enhancer = i[0].split(':')
            enhancer = str(enhancer[0] + ':' + str(enhancer[1]) + ':' + str(enhancer[2]))
            enhancer = enhancer.strip('chr')
            all_ungapped = int(i[4])
            all_identical = int(i[5])
            all_length = int(i[8])
            exclude_ungapped = int(i[6])
            exclude_identical = int(i[7])
            all_exclude = int(i[9])

            align[enhancer] = np.array([all_ungapped, all_identical, all_length, exclude_ungapped,
                                        exclude_identical, all_exclude])

    # Adding non lifted enhancers
    all_enh = {}
    #with open(path_data + "CAGE/CAGE_enhancer_genomic_positions_mm10.bed") as f2:
    with open(path_data + "potential_enhancers/"+enhancers + "_enhancer_genomic_positions_hg38.bed") as f2:
        for i in f2.readlines():
            i = i.strip("\n")
            i = i.split("\t")
            enh = str(i[3]).strip('chr')
            all_enh[enh] = 1        # usefull to get same enh order within species
            if enh not in align.keys():
                align[enh] = np.array([0, 0, 0, 0, 0, 0])

    # Linking enhancer to contacted regions
    enh2frag = {}
    with open(path_data+"overlap/"+origin_sp+"_merged_overlap_" + enhancers + ".txt") as f1:
        for i in f1.readlines()[1:]:
            i = i.strip("\n")
            i = i.split("\t")
            frag = i[0].strip(':+')
            for enhancer in i[4].split(','):
                enhancer = enhancer.strip('chr')
                if enhancer not in enh2frag.keys():
                    enh2frag[enhancer] = [frag]
                else:
                    enh2frag[enhancer].append(frag)

    # Contacted region to baits
    if data == "_simul":
        infile = "/Simulations/simulations_"+origin_sp+"_10Mb_bin5kb_fragoverbin_chr_merged.txt"
    else:
        infile = "/all_interactions/all_interactions_chr_merged.txt_cell_names"

    inter = {}
    with open(path_data+infile) as f3:
        for i in f3.readlines()[1:]:
            i = i.strip("\n")
            i = i.split("\t")
            bait = (str(i[0]) + ":" + str(i[1]) + ":" + str(i[2]))
            merged_PIR = (str(i[3]) + ":" + str(i[4]) + ":" + str(i[5]))
            if i[0] == i[3]:
                dist_obs = float(i[7])
                if 25000 <= dist_obs <= 10000000:
                    if merged_PIR not in inter.keys():
                        inter[merged_PIR] = [(bait, dist_obs)]
                    else:
                        inter[merged_PIR].append((bait, dist_obs))

    # Baits composition
    TSS_count = {}
    genes_count = {}
    with open(path_result + "bait_composition_" + origin_sp + data + "_merged.txt") as f3:
        for i in f3.readlines()[1:]:
            i = i.strip("\n")
            i = i.split("\t")
            bait = (str(i[0]) + ":" + str(i[1]) + ":" + str(i[2]))
            TSS_count[bait] = int(i[len(i)-2])
            genes_count[bait] = int(i[len(i)-1])

    # Writting output
    output = open(path_result + "Sequence_conservation/"+ origin_sp + "/" + enhancers + "/"+origin_sp+"2"+target_sp+data+"_merged.txt", 'w')
    if os.stat(path_result + "Sequence_conservation/"+ origin_sp + "/" + enhancers + "/"+origin_sp+"2"+target_sp+data+"_merged.txt").st_size == 0:
        output.write("enhancer\tall_ungapped\tall_identical\tall_length\texclude_ungapped\texclude_identical"
                     "\tall_exclude\tnb_bait\tnb_TSS\tnb_genes\n")

    for enhancer in all_enh.keys():
        nb_TSS = nb_genes = nb_bait = 0
        if enhancer in enh2frag.keys():     # enhancer can be not present in contacted fragment
            for frag in enh2frag[enhancer]:
                if frag in inter.keys():   # contacted fragment can be not present in filtered interactions (cis + dist)
                    nb_bait = len(inter[frag])
                    for bait in inter[frag]:
                        nb_TSS += TSS_count[bait]
                        nb_genes += genes_count[bait]

        output.write(str(enhancer) + '\t' + str(align[enhancer][0]) + '\t' + str(align[enhancer][1])
                     + '\t' + str(align[enhancer][2]) + '\t' + str(align[enhancer][3])
                     + '\t' + str(align[enhancer][4]) + '\t' + str(align[enhancer][5])
                     + '\t' + str(nb_bait) + '\t' + str(nb_TSS) + '\t' + str(nb_genes) + '\n')

    output.close()


for sp in target_sps:
    for enh in enhancers:
        for data in datas:
            conserv_enh(origin_sp, sp, enh, data)
