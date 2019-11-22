#!/usr/bin/env python3
# coding=utf-8

import os
import numpy as np

# Conservation mouse interaction in human:
origin_sp = "mouse"
target_sp = "human"
data = "_simul"  # or "" for observed data
data_target = "_simul"


path_data = "/home/laverre/Documents/Regulatory_Landscape/data/"
path_result = "/home/laverre/Documents/Regulatory_Landscape/result/"
seuil = 0.1


print("Origin sp:", origin_sp, "; data:", data , "; to data :", data_target)

# Align score for each fragment in origin sp in target sp
frag_conserv = {}
with open(path_result + "conservation/alignments/" + origin_sp + "/contacted_seq/AlignmentStatistics_Excluding_all_Exons_" + origin_sp + "2" + target_sp + "_merged_bait.txt") as f1:
    for i in f1.readlines()[1:]:
        i = i.strip("\n")
        i = i.split("\t")
        origin_frag = i[0].strip(':+')
        origin_frag = origin_frag.strip(':-')
        frag_align = i[1].strip(':+')
        frag_align = frag_align.strip(':-')

        score_all_ungapped = int(i[4]) / (int(i[8])+1)

        if score_all_ungapped > seuil:
            frag_conserv[origin_frag] = (frag_align, score_all_ungapped)

print("Score align : done ! ")

# Duplication score
frag_dupli = {}
with open(path_result + "BLAT_duplication/" + origin_sp + "_merged_fragments_duplication_stats.txt") as f1:
    for i in f1.readlines()[1:]:
        i = i.strip("\n")
        i = i.split("\t")
        frag = i[0]
        frag_dupli[frag] = int(i[3]) - 1

print("Score dupli : done ! ")

# Overlap homologous fragment from origin sp to fragment in target sp
overlap_target = {}
with open(path_data + target_sp + "/overlap/" +origin_sp+"2"+target_sp+"_merged_overlap_"+target_sp+"_merged.txt") as f2:
    for i in f2.readlines()[1:]:
        i = i.strip("\n")
        i = i.split("\t")
        frag_align = str(i[1]+':' + str(i[2]) + ':' + str(i[3]))
        overlap_target[frag_align] = i[4].split(',')  # Add all overlapping fragment

print("Overlap aligned with target sp frag : done !")

target_interaction = {}
stats_target = {}
# Interaction in target sp
if data_target == "_simul":
    input = "/Simulations/simulations_" + target_sp + "_10Mb_bin5kb_fragoverbin_chr_merged.txt"
else:
    input = "/all_interactions/all_interactions_chr_merged.txt_cell_names"

with open(path_data+target_sp+input) as f3:
    for i in f3.readlines()[1:]:
        i = i.strip("\n")
        i = i.split("\t")
        bait = (str(i[0]) + ":" + str(i[1]) + ":" + str(i[2]))
        PIR = (str(i[3]) + ":" + str(i[4]) + ":" + str(i[5]))

        if data_target == "":
            nb_type = float(i[9])
            median_strength = np.median(list(map(float, i[10].split(","))))
        else:
            nb_type = median_strength = "NA"

        if i[0] == i[3]:
            target_dist = float(i[7])
            stats_target[bait + '-' + PIR] = (target_dist, nb_type, median_strength)

            if bait in target_interaction.keys():
                target_interaction[bait].append(PIR)
            else:
                target_interaction[bait] = [PIR]

#print("target", len(target_interaction.keys()))
print("Interaction in target sp : done !")
print("Running conservation of interactions... ")

# Interaction in origin sp to interaction in target sp
if data == "_simul":
    input = "/Simulations/simulations_" + origin_sp + "_10Mb_bin5kb_fragoverbin_chr_merged.txt"
else:
    input = "/all_interactions/all_interactions_chr_merged.txt_cell_names"

conserv_inter = {}
stats_origin = {}
bait_length = {}
PIR_length = {}
PIR_list = []
bait_list = []
count = {}
with open(path_data+origin_sp+input) as f3:
    for i in f3.readlines()[1:]:
        i = i.strip("\n")
        i = i.split("\t")
        bait = (str(i[0]) + ":" + str(i[1]) + ":" + str(i[2]))
        PIR = (str(i[3]) + ":" + str(i[4]) + ":" + str(i[5]))
        bait_length[bait] = int(i[2]) - int(i[1])
        PIR_length[PIR] = int(i[5]) - int(i[4])

        if i[7] != "NA":
            target_dist = float(i[7])

        if data == "":
            nb_type = float(i[9])
            median_strength = np.median(list(map(float, i[10].split(","))))
        else:
            nb_type = median_strength = "NA"

        PIR_list.append(PIR)
        bait_list.append(bait)

        if bait not in frag_dupli.keys():
            frag_dupli[bait] = "NA"
        if PIR not in frag_dupli.keys():
            frag_dupli[PIR] = "NA"

        stats_origin[bait + '-' + PIR] = (target_dist, nb_type, median_strength)
        if i[0] == i[3]:
            if 25000 < abs(target_dist) <= 10000000:
                bait_list.append(bait)
                if bait in frag_conserv.keys():
                    if PIR in frag_conserv.keys():
                        for target_bait in overlap_target[frag_conserv[bait][0]]:
                            for target_PIR in overlap_target[frag_conserv[PIR][0]]:
                                if target_bait in target_interaction.keys():
                                    if target_PIR in target_interaction[target_bait]:
                                        conserv_inter[bait + '-' + PIR] = target_bait + '-' + target_PIR
                                    else:
                                        if bait + '-' + PIR not in conserv_inter.keys():
                                            conserv_inter[bait + '-' + PIR] = "conserv_no_interact"
                                else:
                                    count[target_bait] = 1
                                    if bait + '-' + PIR not in conserv_inter.keys():
                                        conserv_inter[bait + '-' + PIR] = "bait_no_interact"
                    else:
                        if bait + '-' + PIR not in conserv_inter.keys():
                            conserv_inter[bait + '-' + PIR] = "PIR_not_conserv"
                else:
                    if bait + '-' + PIR not in conserv_inter.keys():
                        conserv_inter[bait + '-' + PIR] = "bait_not_conserv"
            else:
                if bait + '-' + PIR not in conserv_inter.keys():
                    conserv_inter[bait + '-' + PIR] = "interactions_not_in_range"
        else:
            if bait + '-' + PIR not in conserv_inter.keys():
                conserv_inter[bait + '-' + PIR] = "trans_interactions"


print("origin", len(set(bait_list)))
print("nb bait no interact", len(count.keys()))

for PIR in set(PIR_list):
    if PIR not in frag_conserv.keys():
        frag_conserv[PIR] = ("NA", 0)
for bait in set(bait_list):
    if bait not in frag_conserv.keys():
        frag_conserv[bait] = ("NA", 0)


print("Writting output...")

output = open(path_result + "/conservation/interaction_conservation/"+origin_sp+data+"_to_"+target_sp+data_target+"_"+str(seuil)+"_merged.txt", 'w')
if os.stat(path_result + "/conservation/interaction_conservation/"+origin_sp+data+"_to_"+target_sp+data_target+"_"+str(seuil)+"_merged.txt").st_size == 0:
    output.write("origin_interaction\torigin_dist\torigin_nb_tissu\torigin_strength\ttarget_interaction\t"
                 "target_dist\ttarget_nb_tissu\ttarget_strength\tbait_score\tbait_length\tbait_dupli\t"
                 "PIR_score\tPIR_length\tPIR_dupli\n")

for inter in conserv_inter.keys():
    bait = inter.split('-')[0]
    PIR = inter.split('-')[1]
    if conserv_inter[inter] in stats_target.keys():
        output.write(str(inter) + '\t' + str(stats_origin[inter][0]) + '\t' + str(stats_origin[inter][1]) + '\t' +
                     str(stats_origin[inter][2]) + '\t' + str(conserv_inter[inter]) + '\t' +
                     str(stats_target[conserv_inter[inter]][0]) + '\t' + str(stats_target[conserv_inter[inter]][1])
                     + '\t' + str(stats_target[conserv_inter[inter]][2]) + '\t' + str(frag_conserv[bait][1]) + '\t' +
                     str(bait_length[bait]) + '\t' + str(frag_dupli[bait]) + '\t' + str(frag_conserv[PIR][1]) + '\t' +
                     str(PIR_length[PIR]) + '\t' + str(frag_dupli[PIR]) + '\n')

    else:
        output.write(str(inter) + '\t' + str(stats_origin[inter][0]) + '\t' + str(stats_origin[inter][1]) + '\t' +
                     str(stats_origin[inter][2]) + '\t' + str(conserv_inter[inter]) + '\t' + "NA" + '\t' + "NA" +
                     '\t' + "NA" + '\t' + str(frag_conserv[bait][1]) + '\t' + str(bait_length[bait]) + '\t' +
                     str(frag_dupli[bait]) + '\t' + str(frag_conserv[PIR][1]) + '\t' + str(PIR_length[PIR]) + '\t' +
                     str(frag_dupli[PIR]) + '\n')


print("All done ! ")
output.close()
