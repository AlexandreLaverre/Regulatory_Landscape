#!/usr/bin/env python3
# coding=utf-8

import os

frag_conserv = {}
with open("../../result/alignments/human2mouse/AlignmentStatistics_TBA_human2mouse_withoutnull.txt") as f2:
    for i in f2.readlines()[1:]:
        i = i.strip("\n")
        i = i.split("\t")
        frag = i[0].strip(':+')
        align = i[1].split(':')
        frag_align = (str(align[0]) + ":" + str(align[1]) + ":" + str(align[2]))
        score = int(i[5])/int(i[2])
        midfrag = ((int(align[1]) + int(align[2])) / 2)

        frag_conserv[frag] = (frag_align, midfrag, score)



inter = {}
comp_dist = {}
score = {}
pairs = pairs_conserv = conserv_intra = 0

#Simulations/simulations_human_10Mb_bin5kb_fragoverbin_chr.txt
with open("../../data/human/all_interactions/all_interactions_chr.txt") as f3:
    for i in f3.readlines()[1:]:
        i = i.strip("\n")
        i = i.split("\t")
        midbait = ((int(i[1]) + int(i[2])) / 2)
        midcontact = ((int(i[4]) + int(i[5])) / 2)
        bait = (str(i[0]) + ":" + str(int(i[1])-1) + ":" + str(i[2]))
        PIR = (str(i[3]) + ":" + str(int(i[4])-1) + ":" + str(i[5]))
        dist_obs = abs(midbait - midcontact)

        if i[0] == i[3]:
            if 25000 < abs(midbait - midcontact) <= 10000000:
                pairs += 1
                if bait in frag_conserv.keys():
                    if PIR in frag_conserv.keys():
                        pairs_conserv += 1

                        chr_bait_conserv = str(frag_conserv[bait][0].split(':')[0])
                        chr_PIR_conserv = str(frag_conserv[PIR][0].split(':')[0])

                        if chr_bait_conserv == chr_PIR_conserv:
                            inter[bait+'\t'+PIR] = str(dist_obs) + '\t' + str(frag_conserv[bait][0]) + '\t' + \
                                                   str(frag_conserv[PIR][0]) + '\t' + str(abs(frag_conserv[bait][1] - frag_conserv[PIR][1]))\
                                                   + '\t' + str(frag_conserv[bait][2]) + '\t' + str(frag_conserv[PIR][2])

                            conserv_intra += 1

                        else:
                            inter[bait + '\t' + PIR] = str(dist_obs) + '\t' + str(frag_conserv[bait][0]) + '\t' + \
                                                       str(frag_conserv[PIR][0]) + '\t' + "NA" + '\t' + \
                                                       str(frag_conserv[bait][2]) + '\t' + str(frag_conserv[PIR][2])
                    else:
                        inter[bait + '\t' + PIR] = str(dist_obs) + '\t' + str(frag_conserv[bait][0]) + '\t' + \
                                                   "NA" + '\t' + "NA" + '\t' + str(frag_conserv[bait][2]) + '\t' + "NA"

                else:
                    inter[bait + '\t' + PIR] = str(dist_obs) + '\t' + "NA" + '\t' + "NA" + '\t' + "NA" + '\t' + "NA" + '\t' + "NA"

print("Paires conservées :", pairs_conserv, 'dont_intra', conserv_intra, ' sur ', pairs)


output = open("../../result/alignments/human2mouse/interaction_conserv.txt", 'w')
if os.stat("../../result/alignments/human2mouse/interaction_conserv.txt").st_size == 0:
    output.write("bait\tPIR\tdist_obs\tbait_conserv\tPIR_conserv\tdist_conserv\tbait_score\tPIR_score\n")

for k, v in inter.items():
    output.write(k + '\t' + v + '\n')

output.close()

