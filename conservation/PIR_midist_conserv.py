#!/usr/bin/env python3
# coding=utf-8

import os
import numpy as np

frag_conserv = {}
with open("../../result/alignments/mouse2human/AlignmentStatistics_TBA_mouse2human_withoutnull.txt") as f2:
    for i in f2.readlines()[1:]:
        i = i.strip("\n")
        i = i.split("\t")
        frag = i[0].split(':')
        frag = str(frag[0]+':' + str(int(frag[1])+1) + ':' + str(frag[2]))
        align = i[1].split(':')
        frag_align = (str(align[0]) + ":" + str(int(align[1])+1) + ":" + str(align[2]))
        score = int(i[5])/int(i[2])
        midfrag = ((int(align[1])+1 + int(align[2])) / 2)

        frag_conserv[frag] = (frag_align, midfrag, score)

exon_pb = {}
with open("../../data/mouse/overlap/mouse_frag_overlap_exons.txt") as f2:
    for i in f2.readlines()[1:]:
        i = i.strip("\n")
        i = i.split("\t")
        frag = ('chr' + str(i[0]) + ':' + str(i[1]) + ':' + str(i[2]))

        length = str(int(i[4]))
        pb = str(int(i[5]))

        exon_pb[frag] = str(length + '\t' + pb)

#
#Simulations/simulations_mouse_10Mb_bin5kb_fragoverbin_chr.txt
inter = {}
conserv = {}
with open("../../data/mouse/all_interactions/all_interactions_chr.txt") as f3:
    for i in f3.readlines()[1:]:
        i = i.strip("\n")
        i = i.split("\t")
        midbait = ((int(i[1]) + int(i[2])) / 2)
        midcontact = ((int(i[4]) + int(i[5])) / 2)
        bait = (str(i[0]) + ":" + str(i[1]) + ":" + str(i[2]))
        PIR = (str(i[3]) + ":" + str(i[4]) + ":" + str(i[5]))
        dist_obs = abs(midbait - midcontact)

        if i[0] == i[3]:
            if 25000 < abs(midbait - midcontact) <= 10000000:

                if PIR in inter.keys():
                    inter[PIR].append(dist_obs)
                else:
                    inter[PIR] = [dist_obs]

                if PIR in frag_conserv.keys():
                    conserv[PIR] = str(frag_conserv[PIR][2])
                else:
                    conserv[PIR] = str(0)


output = open("../../result/alignments/mouse2human/PIR_midist_conserv.txt", 'w')
if os.stat("../../result/alignments/mouse2human/PIR_midist_conserv.txt").st_size == 0:
    output.write("PIR\tnb_contact\tmidist_obs\tPIR_score\tlength\texon_pb\n")

for PIR in inter.keys():
    output.write(PIR + '\t' + str(len(inter[PIR])) + '\t' + str(np.median(inter[PIR])) + '\t' + conserv[PIR]
                 + '\t' + exon_pb[PIR] + '\n')

output.close()

