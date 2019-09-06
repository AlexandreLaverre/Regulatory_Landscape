#!/usr/bin/env python3
# coding=utf-8

import os
import sys
import re

path_data = "/home/laverre/Documents/Regulatory_Landscape/data/mouse/Simulations/"

reference_file = path_data + "PIR.txt"
interest_file = path_data + "bait.txt"
output_file = path_data + "mouse_PIR_overlap_bait_simul.txt"


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
            pos = (int(i[1]), int(i[2]), str(i[0]))  # value = position (start + end)

            if str(i[0]) in dic.keys():
                dic[str(i[0])].append(pos)     # keys = chromosome
            else:
                dic[str(i[0])] = [pos]

    # Sorting tuples for each key
    for k in dic.keys():
        dic[k] = list(set(tuple(x) for x in dic[k]))
        dic[k].sort(key=lambda x: x[0])

    return dic


ref_dic = sorted_dictionary(reference_file)
int_dic = sorted_dictionary(interest_file)
print("Dictionaries ready")

print("Running overlap... ")

# Overlap interest to reference
dic_output = {}
for chr in ref_dic.keys():
    first_i = 0
    for pos in ref_dic[chr]:
        start = pos[0]
        end = pos[1]
        ref_pos = str(chr) + "\t" + str(start) + "\t" + str(end) #+ "\t" + str(pos[2])

        if chr in int_dic.keys():

            # Initialization of first possible overlapping interest position
            i = first_i
            while i < len(int_dic[chr]) and int_dic[chr][i][1] < start:  # -250:
                i += 1
            first_i = i

            # Adding all overlapping interest position to reference position
            while i < len(int_dic[chr]) and int_dic[chr][i][0] <= end:  #+ 250:
                if ref_pos in dic_output.keys():
                    if any(int_dic[chr][i][2] in overlap for overlap in dic_output[ref_pos]):
                        a = 1
                    else:
                        dic_output[ref_pos].append(int_dic[chr][i])
                else:
                    dic_output[ref_pos] = [int_dic[chr][i]]
                i += 1

            # Adding reference position without overlap
            if ref_pos not in dic_output.keys():
                dic_output[ref_pos] = [('NA', 'NA', 'NA')]

        else:
            dic_output[ref_pos] = [('NA', 'NA', 'NA')]

print("Writting output... ")

output = open(output_file, 'w')
if os.stat(output_file).st_size == 0:
    output.write("chr\tstart\tend\toverlap_ID\n")

for ref_pos, int_pos in dic_output.items():
    output.write(ref_pos + "\t")   # add chr if not chr in ID ref
    chr = ref_pos.split('\t')[0]
    count = 0

    for i in int_pos:
        count += 1
        if count == len(int_pos):
            if str(i[0]) == 'NA':
                output.write('unbait' + "\n")
            else:
                output.write(str(i[2]) + ':' + str(i[0]) + ':' + str(i[1]) + "\n")
        else:
            output.write(str(i[2]) + ':' + str(i[0]) + ':' + str(i[1]) + ",")  # Only ID

output.close()
print("All done !")
