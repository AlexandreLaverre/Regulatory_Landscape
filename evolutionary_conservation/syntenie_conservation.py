#!/usr/bin/env python3
# coding=utf-8

import os
import numpy as np


def sorted_dictionary(origin_sp, target_sp, data):
    print(origin_sp, "to", target_sp, ":", data)
    # Align score for each fragment in origin sp in target sp
    frag_conserv = {}
    with open("../../result/conservation/alignments/"+origin_sp+"/contacted_seq/AlignmentStatistics_Excluding_all_Exons_"+origin_sp+"2"+target_sp+"_merged_bait.txt") as f1:
        for i in f1.readlines()[1:]:
            i = i.strip("\n")
            i = i.split("\t")
            origin_frag = i[0].strip(':+')
            origin_frag = origin_frag.strip(':-')
            frag_align = i[1].strip(':+')
            frag_align = frag_align.strip(':-')

            score = int(i[6]) / (int(i[9]) + 1)  # exclude_ungapped

            frag_conserv[origin_frag] = (frag_align, score)

    # Duplication score
    frag_dupli = {}
    with open("../../result/BLAT_duplication/" + origin_sp + "_merged_fragments_duplication_stats.txt") as f1:
        for i in f1.readlines()[1:]:
            i = i.strip("\n")
            i = i.split("\t")
            frag = i[0]
            frag_dupli[frag] = int(i[3]) - 1

    output = open("../../result/conservation/syntenie_conservation/"+origin_sp+"/"+origin_sp+"2"+target_sp+"_conservation_syntenie_with_notconserv"+data+".txt", 'w')
    if os.stat("../../result/conservation/syntenie_conservation/"+origin_sp+"/"+origin_sp+"2"+target_sp+"_conservation_syntenie_with_notconserv"+data+".txt").st_size == 0:
        output.write("origin_interaction\torigin_dist\tnb_type\tstrength\tbait_lift\tbait_score\tbait_dupli\tbait_length\t"
                     "PIR_lift\tPIR_score\tPIR_dupli\tPIR_length\ttarget_dist\n")

    # Interaction in origin sp
    if data == "_simul":
        infile = "/Simulations/simulations_"+origin_sp+"_10Mb_bin5kb_fragoverbin_chr_merged.txt"
    else:
        infile = "/all_interactions/all_interactions_chr_merged.txt_cell_names"

    dic_bait = {}
    dic_PIR = {}
    dic_bait_dupli = {}
    dic_bait_nodupli = {}
    dic_PIR_dupli = {}
    dic_PIR_nodupli = {}

    with open("../../data/"+origin_sp+infile) as f3:
        for i in f3.readlines()[1:]:
            i = i.strip("\n")
            i = i.split("\t")
            midbait = ((int(i[1]) + int(i[2])) / 2)
            bait_length = int(i[2]) - int(i[1])

            midcontact = ((int(i[4]) + int(i[5])) / 2)
            PIR_length = int(i[5]) - int(i[4])

            bait = (str(i[0]) + ":" + str(i[1]) + ":" + str(i[2]))
            PIR = (str(i[3]) + ":" + str(i[4]) + ":" + str(i[5]))
            origin_dist = abs(midbait - midcontact)

            median_strength = i[10]
            nb_type = i[9]

            if i[0] == i[3]:
                if 25000 < abs(midbait - midcontact) <= 10000000:
                    if bait in frag_dupli.keys():
                        bait_dupli = frag_dupli[bait]
                        dic_bait_dupli[bait] = 1
                    else:
                        bait_dupli = "NA"
                        dic_bait_nodupli[bait] = 1

                    dic_bait[bait] = 1
                    dic_PIR[PIR] = 1

                    if PIR in frag_dupli.keys():
                        PIR_dupli = frag_dupli[PIR]
                        dic_PIR_dupli[PIR] = 1
                    else:
                        PIR_dupli = "NA"
                        dic_PIR_nodupli[PIR] = 1

                    if bait in frag_conserv.keys():
                        bait_lift = frag_conserv[bait]
                        if PIR in frag_conserv.keys():
                            PIR_lift = frag_conserv[PIR]

                            if bait_lift[0].split(':')[0] == PIR_lift[0].split(':')[0]:
                                midbait_lift = ((int(bait_lift[0].split(':')[1]) + int(bait_lift[0].split(':')[2])) / 2)
                                midPIR_lift = ((int(PIR_lift[0].split(':')[1]) + int(PIR_lift[0].split(':')[2])) / 2)
                                target_dist = abs(midbait_lift - midPIR_lift)

                            else:
                                target_dist = "NA"

                            output.write(bait + '-' + PIR + '\t' + str(origin_dist) + '\t' + str(nb_type) + '\t' +
                                         str(median_strength) + '\t' + bait_lift[0] + '\t' + str(bait_lift[1]) + '\t' +
                                         str(bait_dupli) + '\t' + str(bait_length) + '\t' + PIR_lift[0] + '\t' +
                                         str(PIR_lift[1]) + '\t' + str(PIR_dupli) + '\t' + str(PIR_length) + '\t' +
                                         str(target_dist) + '\n')

                        else:
                            output.write(bait + '-' + PIR + '\t' + str(origin_dist) + '\t' + str(nb_type) + '\t' +
                                         str(median_strength) + '\t' + bait_lift[0] + '\t' + str(bait_lift[1]) + '\t' +
                                         str(bait_dupli) + '\t' + str(bait_length) + '\t' + "NA" + '\t' + "NA" +
                                         '\t' + str(PIR_dupli) + '\t' + str(PIR_length) + '\t' + "NA" + '\n')
                    else:
                        output.write(bait + '-' + PIR + '\t' + str(origin_dist) + '\t' + str(nb_type) + '\t' +
                                     str(median_strength) + '\t' + "NA" + '\t' + "NA" + '\t' + str(bait_dupli) + '\t' +
                                     str(bait_length) + '\t' + "NA" + '\t' + "NA" + '\t' + str(PIR_dupli) + '\t' +
                                     str(PIR_length) + '\t' + "NA" + '\n')

        print("nb_bait", len(dic_bait.keys()))
        print("nb_PIR", len(dic_PIR.keys()))
        print("nb_bait_dupli", len(dic_bait_dupli.keys()))
        print("nb_bait_nodupli", len(dic_bait_nodupli.keys()))
        print("nb_PIR_dupli", len(dic_PIR_dupli.keys()))
        print("nb_PIR_nodupli", len(dic_PIR_nodupli.keys()))

    output.close()


origin_sp = "human"
target_sps = ["chicken"] #, "opossum", "rat", "mouse", "rabbit", "dog", "elephant", "cow", "macaque"]
datas = ["_simul", ""]

for sp in target_sps:
    for data in datas:
        sorted_dictionary(origin_sp, sp, data)


print("All done ! ")
