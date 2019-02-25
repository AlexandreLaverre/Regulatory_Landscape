#!/usr/bin/env python3
# coding=utf-8

import os

bait_conserv = {}
with open("../../result/alignments/mouse2human/Align_baited_fragments_withoutnull.txt") as f1:
    for i in f1.readlines():
        i = i.split("\t")
        bait = i[0].strip(':+')
        align = i[1].split(':')
        bait_align = (str(align[0]) + ":" + str(align[1]) + ":" + str(align[2]))
        score = int(i[5])/int(i[2])
        midbait = ((int(align[1]) + int(align[2])) / 2)

        bait_conserv[bait] = (bait_align, midbait, score)


PIR_conserv = {}
with open("../../result/alignments/mouse2human/Align_contacted_fragments_withoutnull.txt") as f2:
    for i in f2.readlines():
        i = i.split("\t")
        PIR = i[0].strip(':+')
        align = i[1].split(':')
        PIR_align = (str(align[0]) + ":" + str(align[1]) + ":" + str(align[2]))
        score = int(i[5])/int(i[2])
        midPIR = ((int(align[1]) + int(align[2])) / 2)

        PIR_conserv[PIR] = (PIR_align, midPIR, score)


inter = {}
comp_dist = {}
score = {}
pairs = pairs_conserv = conserv_intra = 0

with open("../../data/mouse/all_interactions/all_interactions_chr.txt") as f3:
    for i in f3.readlines()[1:]:
        i = i.split("\t")
        midbait = ((int(i[1]) + int(i[2])) / 2)
        midcontact = ((int(i[4]) + int(i[5])) / 2)
        bait = (str(i[0]) + ":" + str(int(i[1])-1) + ":" + str(i[2]))
        PIR = (str(i[3]) + ":" + str(int(i[4])-1) + ":" + str(i[5]))
        dist_obs = abs(midbait - midcontact)

        if i[0] == i[3]:
            pairs += 1
            if bait in bait_conserv.keys():
                if PIR in PIR_conserv.keys():
                    pairs_conserv += 1

                    chr_bait_conserv = str(bait_conserv[bait][0].split(':')[0])
                    chr_PIR_conserv = str(PIR_conserv[PIR][0].split(':')[0])

                    if chr_bait_conserv == chr_PIR_conserv:
                        inter[bait+'\t'+PIR] = str(dist_obs) + '\t' + str(bait_conserv[bait][0]) + '\t' + \
                                               str(PIR_conserv[PIR][0]) + '\t' + str(abs(bait_conserv[bait][1] - PIR_conserv[PIR][1]))\
                                               + '\t' + str(bait_conserv[bait][2]) + '\t' + str(PIR_conserv[PIR][2])

                        conserv_intra += 1

                    else:
                        inter[bait + '\t' + PIR] = str(dist_obs) + '\t' + str(bait_conserv[bait][0]) + '\t' + \
                                                   str(PIR_conserv[PIR][0]) + '\t' + "NA" + '\t' + \
                                                   str(bait_conserv[bait][2]) + '\t' + str(PIR_conserv[PIR][2])


print("Paires conservées :", pairs_conserv, 'dont_intra', conserv_intra, ' sur ', pairs)


output = open("../../result/alignments/mouse2human/interaction_conserv.txt", 'w')
if os.stat("../../result/alignments/mouse2human/interaction_conserv.txt").st_size == 0:
    output.write("bait\tPIR\tdist_obs\tbait_conserv\tPIR_conserv\tdist_conserv\tbait_score\tPIR_score\n")

for k, v in inter.items():
    output.write(k + '\t' + v + '\n')

output.close()

