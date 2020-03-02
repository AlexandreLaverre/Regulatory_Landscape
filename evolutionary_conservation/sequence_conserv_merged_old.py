#!/usr/bin/env python3
# coding=utf-8

import os
import numpy as np

def conserv_seq(origin_sp, target_sp, data):
    print(origin_sp, "to", target_sp, ":", data)

    # Alignment scores
    def dic_score(file):
        frag_conserv = {}
        with open("../../result/conservation/alignments/Alignments_statistics/"+origin_sp+"/contacted_seq/"+file) as f2:
            for i in f2.readlines()[1:]:
                i = i.strip("\n")
                i = i.split("\t")
                frag = i[0].split(':')
                frag = str(frag[0]+':' + str(frag[1]) + ':' + str(frag[2]))
                all_ungapped = int(i[4])
                all_identical = int(i[5])
                all_length = int(i[8])
                exclude_ungapped = int(i[6])
                exclude_identical = int(i[7])
                all_exclude = int(i[9])

                frag_conserv[frag] = np.array([all_ungapped, all_identical, all_length, exclude_ungapped,
                                               exclude_identical, all_exclude])

        return frag_conserv

    all_exon = "AlignmentStatistics_Excluding_all_Exons_"+origin_sp+"2"+target_sp+"_merged.txt"
    score_extract_all = dic_score(all_exon)

    if data == "_simul":
        infile = "/Simulations/simulations_"+origin_sp+"_10Mb_bin5kb_fragoverbin_chr"+merge+".txt"
    else:
        infile = "/all_interactions/all_interactions_chr"+merge+".txt_cell_names"

    inter = {}
    with open("../../data/"+origin_sp+infile) as f3:
        for i in f3.readlines()[1:]:
            i = i.strip("\n")
            i = i.split("\t")
            bait = (str(i[0]) + ":" + str(i[1]) + ":" + str(i[2]))
            merged_PIR = (str(i[3]) + ":" + str(i[4]) + ":" + str(i[5]))

            if i[0] == i[3]:
                dist_obs = float(i[7])
                if 25000 <= dist_obs <= 10000000:

                    if merged_PIR not in inter.keys():
                        inter[merged_PIR] = [dist_obs]
                    else:
                        inter[merged_PIR].append(dist_obs)
                    # Adding non-conserved fragments
                    if bait not in score_extract_all.keys():
                        score_extract_all[bait] = np.array([0, 0, 0, 0, 0, 0])

                    if merged_PIR not in score_extract_all.keys():
                        score_extract_all[merged_PIR] = np.array([0, 0, 0, 0, 0, 0])


    print("Writting output...")
    output = open("../../result/conservation/Sequence_conservation/PIR_cons_all_overlap_PECAN_"+origin_sp+"2"+target_sp+data+merge+".txt_onlyconserv", 'w')
    if os.stat("../../result/conservation/Sequence_conservation/PIR_cons_all_overlap_PECAN_"+origin_sp+"2"+target_sp+data+merge+".txt_onlyconserv").st_size == 0:
        output.write("chr\tstart\tend\tall_ungapped\tall_identical\tall_length\texclude_ungapped\texclude_identical\tall_exclude\n")

    for PIR in inter.keys():
        output.write(PIR.split(':')[0] + '\t' + PIR.split(':')[1] + '\t' + PIR.split(':')[2]
                     + '\t' + str(score_extract_all[PIR][0]) + '\t' + str(score_extract_all[PIR][1])
                     + '\t' + str(score_extract_all[PIR][2]) + '\t' + str(score_extract_all[PIR][3])
                     + '\t' + str(score_extract_all[PIR][4]) + '\t' + str(score_extract_all[PIR][5]) + '\n')

    output.close()


origin_sp = "human"
target_sps = ["chicken", "opossum", "rat", "mouse", "rabbit", "dog", "elephant", "cow", "macaque"]
datas = ["", "_simul"]
merge = "_merged"

for sp in target_sps:
    for data in datas:
        conserv_seq(origin_sp, sp, data)

print("All done ! ")
