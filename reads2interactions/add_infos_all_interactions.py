#!/usr/bin/env python3
# coding=utf-8

import re
import os
import numpy as np

# Conservation mouse interaction in human:
sp = "mouse"
data = ""  # or "_simul"
path_data = "/home/laverre/Documents/Regulatory_Landscape/data/"+sp+"/"

if data == "_simul":
    data_type = "simulated"
    infile = path_data + "/Simulations/simulations_" + sp + "_10Mb_bin5kb_fragoverbin_chr.txt"
    output_file = path_data + "/Simulations/simulations_" + sp + "_10Mb_bin5kb_fragoverbin_chr_info.txt"
else:
    data_type = "observed"
    infile = path_data + "/all_interactions/all_interactions_chr.txt"
    output_file = path_data + "/all_interactions/all_interactions_chr_info.txt"

print("specie:", sp)
print("data:", data_type)


# Interactions dict
def sorted_dictionary(file):
    dic = {}
    with open(file, 'r') as f:
        for i in f.readlines()[1:]:
            i = i.strip("\n")
            i = i.split("\t")
            bait = (str(i[0]) + "\t" + str(i[1]) + "\t" + str(i[2]))

            # Nb cell type
            if data == "":
                contact_strength = [float(x) for x in i[6:] if x != "NA"]
                median_strength = np.median(contact_strength)

                if sp == "mouse":
                    types = [i[6], min(i[7:11]), min(i[11], i[12]), i[13], min(i[14], i[15]), i[16], min(i[17:20])]
                    nb_type = len([float(x) for x in types if x != "NA"])

                elif sp == "human":
                    types = [min(i[6], i[22], i[30]), i[7], min(i[8:11]), i[11], i[12], i[13], i[14], i[15], i[16],
                             i[17], min(i[18], i[31]), i[19], i[20], i[21], min(i[23], i[24], i[25], i[29]),
                             min(i[26], i[27], i[28])]
                    nb_type = len([float(x) for x in types if x != "NA"])

            else:
                median_strength = "NA"
                nb_type = "NA"

            PIR = (i[3], i[4], i[5], 'NA', nb_type, median_strength)
            if bait not in dic.keys():
                dic[bait] = [PIR]
            else:
                dic[bait].append(PIR)     # keys = bait[contacted_chr]

    return dic


new_dic = sorted_dictionary(infile)

# PIR overlap bait
PIR_bait = {}
with open("../../data/" + sp + "/" + sp + "_PIR_overlap_bait" + data + ".txt") as f1:
    for i in f1.readlines()[1:]:
        i = i.strip("\n")
        i = i.split("\t")
        frag = str(i[0] + ':' + str(i[1]) + ':' + str(i[2]))
        if i[3] == "unbait":
            PIR_bait[frag] = "unbaited"
        else:
            PIR_bait[frag] = "baited"

print("Writting output... ")

output = open(output_file, 'w')
if os.stat(output_file).st_size == 0:
    output.write("bait_chr\tbait_start\tbait_end\tchr\tstart\tend\tbaited_frag\tdist\tmerged_info\tnb_type\tmed_strength\tnb_contact\n")

for bait, contacted in new_dic.items():
    for i in contacted:
        if bait.split('\t')[0] == i[0]:
            midbait = ((int(bait.split('\t')[1]) + int(bait.split('\t')[2])) / 2)
            midcontact = ((int(i[1]) + int(i[2])) / 2)
            dist_obs = abs(midbait - midcontact)
        else:
            dist_obs = "NA"

        frag = str(i[0] + ':' + str(i[1]) + ':' + str(i[2]))
        if data == "":
            output.write(bait + "\t" + i[0] + "\t" + i[1] + "\t" + i[2] + "\t" + str(PIR_bait[frag]) + "\t" +
                         str(dist_obs) + "\t" + i[3] + "\t" + str(i[4]) + "\t" + str(i[5]) + "\t" + str(len(new_dic[bait])) + "\n")
        else:
            output.write(bait + "\t" + i[0] + "\t" + i[1] + "\t" + i[2] + "\t" + str(PIR_bait[frag]) + "\t" +
                         str(dist_obs) + "\t" + i[3] + "\t" + "NA" + "\t" + "NA" + "\t" + str(len(new_dic[bait])) + "\n")

output.close()
print("All done !")
