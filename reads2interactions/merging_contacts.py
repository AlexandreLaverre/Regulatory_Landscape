#!/usr/bin/env python3
# coding=utf-8

import re
import os
import numpy as np

# Conservation mouse interaction in human:
sp = "human"
data = "_simul_sample"  # or "_simul"
path_data = "/home/laverre/Documents/Regulatory_Landscape/data/" + sp + "/"

if data == "_simul":
    data_type = "simulated"
    infile = path_data + "/Simulations/simulations_" + sp + "_10Mb_bin5kb_fragoverbin_chr.txt"
    output_file = path_data + "/Simulations/simulations_" + sp + "_10Mb_bin5kb_fragoverbin_chr_merged.txt_corrected2"
elif data == "_simul_sample":
    data_type = "sample"
    infile = path_data + "/Simulations/samples_simulated_all_interactions_bait_all.txt"
    output_file = path_data + "/Simulations/samples_simulated_all_interactions_bait_all_merged.txt_corrected2"
else:
    data_type = "observed"
    infile = path_data + "/all_interactions/all_interactions_chr.txt"
    output_file = path_data + "/all_interactions/all_interactions_chr_merged.txt_cell_names_corrected2"


print("specie:", sp)
print("data:", data_type)

# Baited fragment
PIR_bait = {}
with open(path_data + "/Digest_HindIII_None.txt.baitmap") as f1:
    for i in f1.readlines():
        i = i.strip("\n")
        i = i.split("\t")
        frag = "chr" + str(i[0]) + "\t" + str(i[1]) + "\t" + str(i[2])
        PIR_bait[frag] = "bait"


# Interactions dict
def sorted_dictionary(file):
    dic = {}
    with open(file, 'r') as f:
        first_line = f.readline().strip("\n")
        first_line = first_line.split("\t")

        cell_name = first_line[6:len(first_line)]
        baited = 0
        for i in f.readlines():
            i = i.strip("\n")
            i = i.split("\t")
            bait = (str(i[0]) + "\t" + str(i[1]) + "\t" + str(i[2]))
            midbait = (int(i[1]) + int(i[2])) / 2
            midfrag = (int(i[4]) + int(i[5])) / 2
            dist = abs(midbait - midfrag)

            if 25000 < dist < 10000000:
                if data != "_simul":
                    # Nb cell type
                    if sp == "mouse":
                        # cell_value = i[6:len(first_line)]

                        # Cell types
                        cell_value = [i[6], min(i[7:11]), min(i[11], i[12]), i[13], min(i[14], i[15]), i[16], min(i[17:20])]
                        cell_name = ["EpiSC", "ESC", "ESd", "FLC", "preB", "TSC", "preadipo"]

                    elif sp == "human":
                        # cell_value = i[6:len(first_line)]
                        # Cell types
                        cell_value = [min(i[6], i[22], i[30]), i[7], min(i[8:11]), i[11], i[12], i[13], i[14], i[15], i[16],
                                      i[17], min(i[18], i[31]), i[19], i[20], i[21], min(i[23], i[24], i[25], i[29]),
                                      min(i[26], i[27], i[28])]
                        cell_name = ["Bcell", "CD34", "PEK", "pre_adipo", "cardio", "hESC", "MK", "EP", "Mon", "CD8", "Ery",
                                     "Neu", "Foet", "NB", "TCD4", "Mac"]

                    nb_type = []
                    n = 0
                    for x in cell_value:
                        if x != "NA":
                            nb_type.append(cell_name[n])
                        n += 1
                else:
                    nb_type = ["NA"]

                if data == "":
                    contact_strength = [float(x) for x in i[6:] if x != "NA"]
                    median_strength = np.median(contact_strength)
                else:
                    median_strength = "NA"
                    # nb_type = ["NA"]

                contacted = str(i[3]) + "\t" + str(i[4]) + "\t" + str(i[5])
                if contacted in PIR_bait.keys():
                    bait_status = "baited"
                    baited += 1
                else:
                    bait_status = "unbaited"

                PIR = (i[3], i[4], i[5], 'NA', nb_type, median_strength, bait_status)

                if bait not in dic.keys():
                    dic[bait] = {str(i[3]): [PIR]}
                elif str(i[3]) not in dic[bait].keys():
                    dic[bait][str(i[3])] = [PIR]
                else:
                    dic[bait][str(i[3])].append(PIR)  # keys = bait[contacted_chr]

    print("Nb baited :", baited)
    # Sorting tuples for each key
    for bait in dic.keys():
        for contacted_chr in dic[bait]:
            dic[bait][contacted_chr] = list(tuple(x) for x in dic[bait][contacted_chr])
            dic[bait][contacted_chr].sort(key=lambda x: x[1])

    return dic, cell_name


dic, cell_name = sorted_dictionary(infile)

print("Merging adjacent contacted fragments... \n")
# Merging contacted fragments
new_dic = {}
count = 0
for bait in dic.keys():
    new_dic[bait] = []
    for contacted_chr in dic[bait].keys():
        current_start = dic[bait][contacted_chr][0][1]
        current_end = dic[bait][contacted_chr][0][2]
        current_nb_type = dic[bait][contacted_chr][0][4]
        current_all_types = len(current_nb_type)
        current_strength = dic[bait][contacted_chr][0][5]
        current_bait_status = dic[bait][contacted_chr][0][6]
        current_ID = (str(contacted_chr) + ':' + str(current_start) + ':' + str(current_end))

        if len(dic[bait][contacted_chr]) == 1:
            new_dic[bait].append((contacted_chr, current_start, current_end, current_ID, current_nb_type,
                                  current_strength, current_bait_status, current_all_types))
        else:

            for i in range(1, len(dic[bait][contacted_chr])):
                new_start = dic[bait][contacted_chr][i][1]
                new_end = dic[bait][contacted_chr][i][2]
                new_nb_type = dic[bait][contacted_chr][i][4]
                new_all_types = len(new_nb_type)
                new_strength = dic[bait][contacted_chr][i][5]
                new_bait_status = dic[bait][contacted_chr][i][6]
                new_ID = (str(contacted_chr) + ':' + str(new_start) + ':' + str(new_end))

                if int(new_start) == int(current_end) + 1:
                    current_nb_type = list(set(current_nb_type + new_nb_type))
                    current_strength = str(current_strength) + ',' + str(new_strength)
                    current_all_types = str(current_all_types) + ',' + str(new_all_types)
                    current_ID = str(current_ID) + ',' + str(new_ID)
                    if current_bait_status != new_bait_status:
                        current_bait_status = "baited"
                    current_end = new_end
                else:
                    new_dic[bait].append((contacted_chr, current_start, current_end, current_ID, current_nb_type,
                                          current_strength, current_bait_status, current_all_types))
                    current_start = new_start
                    current_end = new_end
                    current_strength = new_strength
                    current_ID = new_ID
                    current_nb_type = new_nb_type
                    current_bait_status = new_bait_status
                    current_all_types = new_all_types

                if i == len(dic[bait][contacted_chr]) - 1:
                    new_dic[bait].append((contacted_chr, current_start, current_end, current_ID, current_nb_type,
                                          current_strength, current_bait_status, current_all_types))

nb_contact = [len(dic[bait][contacted_chr]) for bait in dic.keys() for contacted_chr in dic[bait].keys()]
new_nb_contact = [len(new_dic[bait]) for bait in new_dic.keys()]
nb_elem = sum(nb_contact)
med_contact = np.median(nb_contact)
new_nb_elem = sum(new_nb_contact)
new_med_contact = np.median(new_nb_contact)

print("\t\t\t Before \t After")
print("nb pairs\t", nb_elem, "\t", new_nb_elem)
print("med contact\t", med_contact, "\t\t", new_med_contact)
print("max contact\t", max(nb_contact), "\t\t", max(new_nb_contact))

print("Writting output... ")

output = open(output_file, 'w')
if os.stat(output_file).st_size == 0:
    output.write("bait_chr\tbait_start\tbait_end\tchr\tstart\tend\tbaited_frag\tdist\tmerged_info"
                 "\tnb_type\tmed_strength\tall_nb_types\tcontact_per_bait\tcell_name\t" + '\t'.join(cell_name) + "\n")

for bait, contacted in new_dic.items():
    for i in contacted:
        if bait.split('\t')[0] == i[0]:
            midbait = ((int(bait.split('\t')[1]) + int(bait.split('\t')[2])) / 2)
            midcontact = ((int(i[1]) + int(i[2])) / 2)
            dist_obs = abs(midbait - midcontact)
        else:
            dist_obs = "NA"

        frag = str(i[0] + ':' + str(i[1]) + ':' + str(i[2]))
        # if data == "":
        cell_presence = []
        for cell in cell_name:
            if cell in i[4]:
                cell_presence.append(str(1))
            else:
                cell_presence.append(str(0))

        output.write(bait + "\t" + i[0] + "\t" + i[1] + "\t" + i[2] + "\t" + str(i[6]) + "\t" +
                     str(dist_obs) + "\t" + i[3] + "\t" + str(len(i[4])) + "\t" + str(i[5]) + "\t" + str(i[7]) + "\t" +
                     str(len(new_dic[bait])) + "\t" + str(','.join(i[4])) + "\t" + str('\t'.join(cell_presence)) + "\n")
        # else:
        #    output.write(bait + "\t" + i[0] + "\t" + i[1] + "\t" + i[2] + "\t" + str(i[6]) + "\t" +
        #                 str(dist_obs) + "\t" + i[3] + "\t" + "NA" + "\t" + "NA" + "\t" + str(len(new_dic[bait]))
        #                 + "\t" + "NA" + "\t" + "NA" + "\n")

output.close()
print("All done !")
