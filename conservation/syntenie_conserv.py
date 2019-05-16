#!/usr/bin/env python3
# coding=utf-8

import os
import numpy as np

# Conservation mouse interaction in human:
origin_sp = "human"
target_sp = "mouse"
data = ""  # or "_simul"

# Align score for each fragment in origin sp in target sp
frag_conserv = {}
with open("../../result/alignments/"+origin_sp+"2"+target_sp+"/AlignmentStatistics_TBA_"+origin_sp+"2"+target_sp+"_withoutnull.txt") as f1:
    for i in f1.readlines()[1:]:
        i = i.strip("\n")
        i = i.split("\t")
        origin_frag = i[0].split(':')
        origin_frag = str(origin_frag[0]+':' + str(int(origin_frag[1])+1) + ':' + str(origin_frag[2]))

        align = i[1].split(':')
        frag_align = (str(align[0]) + ":" + str(int(align[1])+1) + ":" + str(align[2]))
        score = int(i[5])/int(i[2])

        frag_conserv[origin_frag] = (frag_align, score)

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


print("Calculating and writting output...")
output = open("../../result/conservation/"+origin_sp+"2"+target_sp+"_conservation_syntenie_with_notconserv"+data+".txt3", 'w')
if os.stat("../../result/conservation/"+origin_sp+"2"+target_sp+"_conservation_syntenie_with_notconserv"+data+".txt3").st_size == 0:
    output.write("origin_interaction\torigin_dist\tnb_type\tstrength\tbait_lift\tbait_score\tbait_dupli\tbait_length\t"
                 "PIR_lift\tPIR_score\tPIR_dupli\tPIR_length\ttarget_dist\n")

# Interaction in origin sp
if data == "_simul":
    input = "/Simulations/simulations_" + origin_sp + "_10Mb_bin5kb_fragoverbin_chr.txt"
else:
    input = "/all_interactions/all_interactions_chr.txt"

with open("../../data/"+origin_sp+input) as f3:
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
            nb_type = "NA"

        if i[0] == i[3]:
            if 25000 < abs(midbait - midcontact) <= 10000000:
                if bait in frag_dupli.keys():
                    bait_dupli = frag_dupli[bait]
                else:
                    bait_dupli = "NA"

                if PIR in frag_dupli.keys():
                    PIR_dupli = frag_dupli[PIR]
                else:
                    PIR_dupli = "NA"

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

output.close()
print("All done ! ")
