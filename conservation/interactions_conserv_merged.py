#!/usr/bin/env python3
# coding=utf-8

import os
import numpy as np

# Conservation mouse interaction in human:
origin_sp = "mouse"
target_sp = "human"
data = ""  # or "_simul"

print("Origin sp:", origin_sp, "; data:", data)

# Align score for each fragment in origin sp in target sp
frag_conserv = {}
with open("../../result/alignments/"+origin_sp+"2"+target_sp+"/AlignmentStatistics_pecan_Excluding_all_Exons.txt") as f1:
    for i in f1.readlines()[1:]:
        i = i.strip("\n")
        i = i.split("\t")
        origin_frag = i[0].split(':')
        origin_frag = str(origin_frag[0]+':' + str(int(origin_frag[1])+1) + ':' + str(origin_frag[2]))

        align = i[1].split(':')
        frag_align = (str(align[0]) + ":" + str(int(align[1])+1) + ":" + str(align[2]))
        score_all_ungapped = int(i[4]) / (int(i[8])+1)

        if score_all_ungapped > 0.4:  # or 0.4
            frag_conserv[origin_frag] = (frag_align, score_all_ungapped)

print("Score align : done ! ")

# Duplication score
frag_dupli = {}
with open("../../data/"+origin_sp+"/"+origin_sp+"_restriction_fragments_duplication_0.8.txt") as f1:
    for i in f1.readlines()[1:]:
        i = i.strip("\n")
        i = i.split("\t")
        frag = i[0].split(':')
        frag = str(frag[0]+':' + str(int(frag[1])+1) + ':' + str(frag[2]))
        frag_dupli[frag] = int(i[3])-1

print("Score dupli : done ! ")

# Overlap homologous fragment from origin sp to fragment in target sp
overlap_target = {}
with open("../../data/"+target_sp+"/overlap/"+origin_sp+"2"+target_sp+"_frag_overlap_merged_"+target_sp+".txt") as f2:
    for i in f2.readlines()[1:]:
        i = i.strip("\n")
        i = i.split("\t")
        frag_align = str(i[0]+':' + str(i[1]) + ':' + str(i[2]))
        overlap_target[frag_align] = i[3].split(',')  # Add all overlapping fragment

print("Overlap aligned with target sp frag : done !")

target_interaction = {}
stats_target = {}
# Interaction in target sp
with open("../../data/"+target_sp+"/all_interactions/all_interactions_chr_merged.txt") as f3:
    for i in f3.readlines()[1:]:
        i = i.strip("\n")
        i = i.split("\t")
        bait_length = int(i[2]) - int(i[1])
        PIR_length = int(i[5]) - int(i[4])

        bait = (str(i[0]) + ":" + str(i[1]) + ":" + str(i[2]))
        PIR = (str(i[3]) + ":" + str(i[4]) + ":" + str(i[5]))

        median_strength = np.median(list(map(float, i[10].split(","))))

        # Nb cell type
        if data == "":
            nb_type = np.median(list(map(float, i[9].split(","))))
        else:
            nb_type = "NA"

        if i[0] == i[3]:
            origin_dist = float(i[7])
            stats_target[bait + '-' + PIR] = (origin_dist, nb_type, median_strength)

            if bait in target_interaction.keys():
                target_interaction[bait].append(PIR)
            else:
                target_interaction[bait] = [PIR]

print("Interaction in target sp : done !")

# Interaction in origin sp to interaction in target sp
if data == "_simul":
    input = "/Simulations/simulations_" + origin_sp + "_10Mb_bin5kb_fragoverbin_chr.txt"
else:
    input = "/all_interactions/all_interactions_chr.txt"

conserv_inter = {}
stats_origin = {}
bait_length = {}
PIR_length = {}
PIR_list = []
bait_list = []
with open("../../data/"+origin_sp+input) as f3:
    for i in f3.readlines()[1:]:
        i = i.strip("\n")
        i = i.split("\t")
        midbait = ((int(i[1]) + int(i[2])) / 2)
        midcontact = ((int(i[4]) + int(i[5])) / 2)

        bait = (str(i[0]) + ":" + str(i[1]) + ":" + str(i[2]))
        PIR = (str(i[3]) + ":" + str(i[4]) + ":" + str(i[5]))
        bait_length[bait] = int(i[2]) - int(i[1])
        PIR_length[PIR] = int(i[5]) - int(i[4])
        dist_obs = midbait - midcontact

        contact_strength = [float(x) for x in i[6:] if x != "NA"]
        median_strength = np.median(contact_strength)

        # Nb cell type
        if origin_sp == "mouse" and data == "":
            types = [i[6], min(i[7:11]), min(i[11], i[12]), i[13], min(i[14], i[15]), i[16], min(i[17:20])]
            nb_type = len([float(x) for x in types if x != "NA"])

        elif origin_sp == "human" and data == "":
            types = [min(i[6], i[22], i[30]), i[7], min(i[8:11]), i[11], i[12], i[13], i[14], i[15], i[16], i[17],
                     min(i[18], i[31]), i[19], i[20], i[21], min(i[23], i[24], i[25], i[29]), min(i[26], i[27], i[28])]
            nb_type = len([float(x) for x in types if x != "NA"])

        else:
            contact_strength = "NA"
            median_strength = "NA"
            nb_type = "NA"

        PIR_list.append(PIR)
        bait_list.append(bait)

        if bait not in frag_dupli.keys():
            frag_dupli[bait] = "NA"

        if PIR not in frag_dupli.keys():
            frag_dupli[PIR] = "NA"

        if i[0] == i[3]:
            if 25000 < abs(midbait - midcontact) <= 10000000:
                stats_origin[bait + '-' + PIR] = (dist_obs, nb_type, median_strength)
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
                                    if bait + '-' + PIR not in conserv_inter.keys():
                                        conserv_inter[bait + '-' + PIR] = "bait_no_interact"
                    else:
                        if bait + '-' + PIR not in conserv_inter.keys():
                            conserv_inter[bait + '-' + PIR] = "PIR_not_conserv"
                else:
                    if bait + '-' + PIR not in conserv_inter.keys():
                        conserv_inter[bait + '-' + PIR] = "bait_not_conserv"

print("Conservation of origin interaction : done !")

for PIR in set(PIR_list):
    if PIR not in frag_conserv.keys():
        frag_conserv[PIR] = ("NA", 0)
for bait in set(bait_list):
    if bait not in frag_conserv.keys():
        frag_conserv[bait] = ("NA", 0)

print("Writting output...")

output = open("../../result/conservation/"+origin_sp+"2"+target_sp+"_conservation_interaction_pecan_0.4_merged"+data+".txt", 'w')
if os.stat("../../result/conservation/"+origin_sp+"2"+target_sp+"_conservation_interaction_pecan_0.4_merged"+data+".txt").st_size == 0:
    output.write("origin_interaction\torigin_dist\torigin_nb_tissu\torigin_strength\ttarget_interaction\t"
                 "target_dist\ttarget_nb_tissu\ttarget_strength\tbait_score\tbait_length\tbait_dupli\t"
                 "PIR_score\tPIR_length\tPIR_dupli\n")

for inter in conserv_inter.keys():
    bait = inter.split('-')[0]
    PIR = inter.split('-')[1]
    if conserv_inter[inter] in stats_target.keys():
        output.write(inter + '\t' + str(stats_origin[inter][0]) + '\t' + str(stats_origin[inter][1]) + '\t' +
                     str(stats_origin[inter][2]) + '\t' + str(conserv_inter[inter]) + '\t' +
                     str(stats_target[conserv_inter[inter]][0]) + '\t' + str(stats_target[conserv_inter[inter]][1])
                     + '\t' + str(stats_target[conserv_inter[inter]][2]) + '\t' + str(frag_conserv[bait][1]) + '\t' +
                     str(bait_length[bait]) + '\t' + str(frag_dupli[bait]) + '\t' + str(frag_conserv[PIR][1]) + '\t' +
                     str(PIR_length[PIR]) + '\t' + str(frag_dupli[PIR]) + '\n')

    else:
        output.write(inter + '\t' + str(stats_origin[inter][0]) + '\t' + str(stats_origin[inter][1]) + '\t' +
                     str(stats_origin[inter][2]) + '\t' + str(conserv_inter[inter]) + '\t' + "NA" + '\t' + "NA" +
                     '\t' + "NA" + '\t' + str(frag_conserv[bait][1]) + '\t' + str(bait_length[bait]) + '\t' +
                     str(frag_dupli[bait]) + '\t' + str(frag_conserv[PIR][1]) + '\t' + str(PIR_length[PIR]) + '\t' +
                     str(frag_dupli[PIR]) + '\n')

print("All done ! ")
output.close()
