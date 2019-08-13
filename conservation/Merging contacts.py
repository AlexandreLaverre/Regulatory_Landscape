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
    infile = path_data + "/Simulations/simulations_" + sp + "_10Mb_bin5kb_fragoverbin_chr.txt"
    output_file = path_data + "/Simulations/simulations_" + sp + "_10Mb_bin5kb_fragoverbin_chr_merged.txt"
else:
    infile = path_data + "/all_interactions/all_interactions_chr.txt"
    output_file = path_data + "/all_interactions/all_interactions_chr_merged.txt"

print("Origin sp:", sp, "; data:", data)


# Interactions dict
def sorted_dictionary(file):
    dic = {}
    with open(file, 'r') as f:
        start = 0
        if bool(re.search('tart', f.readline())) is True:      # Check if header is present
            start = 1

        f.seek(0)
        for i in f.readlines()[start:]:
            i = i.strip("\n")
            i = i.split("\t")
            if i[0] == i[3]:
                bait = (str(i[0]) + ":" + str(i[1]) + ":" + str(i[2]))
                PIR = (i[3], i[4], i[5])
                midbait = ((int(i[1]) + int(i[2])) / 2)
                midcontact = ((int(i[4]) + int(i[5])) / 2)

                if 25000 < abs(midbait - midcontact) <= 10000000:
                    if bait in dic.keys():
                        dic[bait].append(PIR)     # keys = bait
                    else:
                        dic[bait] = [PIR]

    # Sorting tuples for each key
    for k in dic.keys():
        dic[k] = list(set(tuple(x) for x in dic[k]))
        dic[k].sort(key=lambda x: x[1])

    return dic


dic = sorted_dictionary(infile)
print("Sorted interactions dictionary ready !")


# Merging contacted fragments
new_dic = {}
count = 0
for bait in dic.keys():
    chr = dic[bait][0][0]
    current_start = dic[bait][0][1]
    current_end = dic[bait][0][2]
    new_dic[bait] = []
    merge_info = "no_merged"

    if len(dic[bait]) == 1:
        new_dic[bait].append((chr, current_start, current_end, merge_info))
    else:

        for i in range(1, len(dic[bait])):
            new_start = dic[bait][i][1]
            new_end = dic[bait][i][2]

            if int(new_start) == int(current_end)+1:
                current_end = new_end
                merge_info = "merged"
            else:
                new_dic[bait].append((chr, current_start, current_end, merge_info))
                current_start = new_start
                current_end = new_end
                merge_info = "no_merged"

            if i == len(dic[bait])-1:
                new_dic[bait].append((chr, current_start, current_end, merge_info))

nb_contact = [len(dic[x]) for x in dic.keys()]
new_nb_contact = [len(new_dic[x]) for x in new_dic.keys()]

nb_elem = sum(nb_contact)
med_contact = np.median(nb_contact)
new_nb_elem = sum(new_nb_contact)
new_med_contact = np.median(new_nb_contact)

print("Before:", nb_elem, "pairs; After :", new_nb_elem)
print("Before:", med_contact, "med contact; After :", new_med_contact)
print("Before:", max(nb_contact), "max contact; After :", max(new_nb_contact))

print("Merging done !")
print("Writting output... ")



output = open(output_file, 'w')
if os.stat(output_file).st_size == 0:
    output.write("chr\tstart\tend\tchr\tstart\tend\tmerged\n")

for bait, contacted in new_dic.items():
    for i in contacted:
        output.write(bait + "\t" + i[0] + "\t" + i[1] + "\t" + i[2] + "\t" + i[3] + "\n")

output.close()
print("All done !")
