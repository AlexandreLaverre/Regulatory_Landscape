#!/usr/bin/env python3
# coding=utf-8

import os
import numpy as np

sp1 = "mouse"
sp2 = "human"

frag_conserv = {}
with open("../../result/alignments/"+sp1+"2"+sp2+"/AlignmentStatistics_TBA_"+sp1+"2"+sp2+"_withoutnull_0.1.txt") as f2:
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


def dic_pb(file):
    elem_pb = {}
    with open("../../data/"+sp1+"/overlap/"+sp1+file) as f2:
        for i in f2.readlines()[1:]:
            i = i.strip("\n")
            i = i.split("\t")
            frag = ('chr' + str(i[0]) + ':' + str(i[1]) + ':' + str(i[2]))

            length = str(int(i[4]))
            pb = str(int(i[5]))

            if file == "_frag_overlap_all_exons.txt":
                elem_pb[frag] = str(length + '\t' + pb)
            else:
                elem_pb[frag] = str(pb)

    return elem_pb


all = "_frag_overlap_all_exons.txt"
all_exon = dic_pb(all)
coding = "_frag_overlap_coding_exons.txt"
coding_exon = dic_pb(coding)
nocoding = "_frag_overlap_nocoding_exons.txt"
nocoding_exon = dic_pb(nocoding)
repeat = "_frag_overlap_repeatmasker.txt"
repeat_pb = dic_pb(repeat)
phastcons = "_frag_overlap_phastcons_noexonic250.txt"
phastcons_pb = dic_pb(phastcons)


def dic_count(file):
    elem_pb = {}
    with open("../../data/"+sp1+"/overlap/"+sp1+file) as f2:
        for i in f2.readlines()[1:]:
            i = i.strip("\n")
            i = i.split("\t")
            frag = ('chr' + str(i[0]) + ':' + str(i[1]) + ':' + str(i[2]))

            if i[3] == "NA":
                count = 0
            else:
                count = len(i[3].split(","))

            elem_pb[frag] = str(count)

    return elem_pb


CAGE = "_frag_overlap_CAGE.txt"
CAGE_pb = dic_count(CAGE)
#ENCODE = "_frag_overlap_ENCODE.txt"
#ENCODE_pb = dic_count(ENCODE)
#GRO_seq = "_frag_overlap_GRO_seq.txt"
#GRO_seq_pb = dic_count(GRO_seq)
#RoadMap = "_frag_overlap_RoadMap.txt"
#RoadMap_pb = dic_count(RoadMap)


#all_interactions/all_interactions_chr.txt
#Simulations/simulations_"+sp1+"_10Mb_bin5kb_fragoverbin_chr.txt
inter = {}
conserv = {}
with open("../../data/"+sp1+"/all_interactions/all_interactions_chr.txt") as f3:
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


output = open("../../result/alignments/"+sp1+"2"+sp2+"/PIR_cons_all_overlap.txt", 'w')
if os.stat("../../result/alignments/"+sp1+"2"+sp2+"/PIR_cons_all_overlap.txt").st_size == 0:
    output.write("PIR\tnb_contact\tmidist_obs\tPIR_score\tlength\tall_exon_pb\tcoding_exon_pb\tnocoding_exon_pb\t"
                 "repeat_pb\tphastcons_noexonic250\tCAGE_pb\n") #\tENCODE_pb\tGRO_seq_pb\tRoadMap_pb\n")

for PIR in inter.keys():
    output.write(PIR + '\t' + str(len(inter[PIR])) + '\t' + str(np.median(inter[PIR])) + '\t' + conserv[PIR]
                 + '\t' + all_exon[PIR] + '\t' + coding_exon[PIR] + '\t' + nocoding_exon[PIR] + '\t' + repeat_pb[PIR]
                 + '\t' + phastcons_pb[PIR] + '\t' + CAGE_pb[PIR] + '\n')
                 #+ '\t' + ENCODE_pb[PIR] + '\t' + GRO_seq_pb[PIR] + '\t' + RoadMap_pb[PIR]

output.close()
