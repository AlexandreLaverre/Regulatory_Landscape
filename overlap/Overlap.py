#!/usr/bin/env python3
# coding=utf-8

import os
import sys
import re

path_data = "/home/laverre/Documents/Regulatory_Landscape/data/mouse/"

reference_file = path_data + sys.argv[1]
interest_file = path_data + sys.argv[2]
output_file = path_data + sys.argv[3]
collapse = sys.argv[4]


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
            pos = (int(i[1]), int(i[2]), str(i[3]))  # value = position (start + end) + ID

            if i[0] in dic.keys():
                dic[i[0]].append(pos)     # keys = chromosome
            else:
                dic[i[0]] = [pos]

    # Sorting tuples for each key
    for k in dic.keys():
        dic[k] = list(set(tuple(x) for x in dic[k]))
        dic[k].sort(key=lambda x: x[0])

    return dic


ref_dic = sorted_dictionary(reference_file)
print("Reference dictionary ready")

int_dic = sorted_dictionary(interest_file)

# Testing overlap in interest dic
if list(int_dic.values())[0][0][0] != list(int_dic.values())[0][0][1]:
    print("Length of interest seq. > 1 pb ")
    if collapse == "T":
        new_int_dic = {}
        for k in int_dic.keys():
            current_start = int_dic[k][0][0]
            current_end = int_dic[k][0][1]
            current_ID = int_dic[k][0][2]

            new_int_dic[k] = []
            for i in range(1, len(int_dic[k])):
                new_start = int_dic[k][i][0]
                new_end = int_dic[k][i][1]
                new_ID = int_dic[k][i][2]
                if current_end >= new_start >= current_start:
                    current_end = new_end
                    current_ID = str(current_ID) + "_" + str(new_ID)
                else:
                    new_int_dic[k].append((current_start, current_end, current_ID))
                    current_start = new_start
                    current_end = new_end
                    current_ID = new_ID

                if i == len(int_dic[k])-1:
                    new_int_dic[k].append((current_start, current_end, current_ID))

        nb_elem = sum(len(int_dic[x]) for x in int_dic.keys())
        new_nb_elem = sum(len(new_int_dic[x]) for x in new_int_dic.keys())
        if nb_elem != new_nb_elem:
            print("There is overlap in interest file !")
        else:
            print("There is no overlap in interest file !")

        int_dic = new_int_dic

else:
    print("Length of interest seq. = 1 pb")


print("Interest dictionary ready")
print("Running overlap... ")

# Overlap interest to reference
dic_output = {}
for chr in ref_dic.keys():
    first_i = 0
    for pos in ref_dic[chr]:
        start = pos[0]
        end = pos[1]
        ref_pos = str(chr) + "\t" + str(start) + "\t" + str(end)

        if chr in int_dic.keys():

            # Initialization of first possible overlapping interest position
            i = first_i
            while i < len(int_dic[chr]) and int_dic[chr][i][1] < start:
                i += 1
            first_i = i

            # Adding all overlapping interest position to reference position
            while i < len(int_dic[chr]) and int_dic[chr][i][0] <= end:
                if ref_pos in dic_output.keys():
                    dic_output[ref_pos].append(int_dic[chr][i])
                else:
                    dic_output[ref_pos] = [int_dic[chr][i]]
                i += 1

            # Adding reference position without overlap
            if ref_pos not in dic_output.keys():
                dic_output[ref_pos] = [('NA', 'NA', 'NA')]

        else:
            dic_output[ref_pos] = [('NA', 'NA', 'NA')]

print("Counting base pair...")
count_bp = {}
length_pos = {}
for ref_pos in dic_output.keys():
    for overlap in dic_output[ref_pos]:
        frag_pos = ref_pos.split('\t')
        length_frag = int(frag_pos[2]) - int(frag_pos[1])
        length_pos[ref_pos] = length_frag
        if overlap[2] != 'NA':
            overlap_start = int(overlap[0])
            overlap_end = int(overlap[1])

            if overlap_start < int(frag_pos[1]):
                overlap_start = int(frag_pos[1])
            if overlap_end > int(frag_pos[2]):
                overlap_end = int(frag_pos[2])

            length_overlap = overlap_end - overlap_start

            if ref_pos in count_bp.keys():
                count_bp[ref_pos] = int(count_bp[ref_pos]) + length_overlap
            else:
                count_bp[ref_pos] = length_overlap

        else:
            count_bp[ref_pos] = 0


print("Writting output... ")

output = open(output_file, 'w')
if os.stat(output_file).st_size == 0:
    output.write("chr\tstart\tend\toverlap_ID\tlength_frag\tnb_bp_overlap\t%\n")

for ref_pos, int_pos in dic_output.items():
    output.write(ref_pos + "\t")
    count = 0

    for i in int_pos:
        count += 1
        if count == len(int_pos):
            output.write(str(i[2]) + "\t" + str(length_pos[ref_pos]) + "\t" + str(count_bp[ref_pos]) + '\t' +
                         str((count_bp[ref_pos]/length_pos[ref_pos])*100) + "\n")
            #output.write(str(i[2])+":"+str(i[0])+"-"+str(i[1])+"\n")
        else:
            output.write(str(i[2]) + ",")
            #output.write(str(i[2])+":"+str(i[0])+"-"+str(i[1]) + ",")

output.close()
print("All done !")
