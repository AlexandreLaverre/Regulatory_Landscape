#!/usr/bin/env python3
# coding=utf-8

import os
import sys


reference_file = sys.argv[1]
interest_file = sys.argv[2]
interest_name = sys.argv[3]
output_file = sys.argv[4]


def sorting_dictionary(file):
    dic = {}
    with open(file, 'r') as f:
        for i in f.readlines():
            i = i.split("\t")
            pos = (int(i[1]), int(i[2]))  # values = positions

            if int(i[1]) == int(i[2]):

            if chr in dic.keys():
                dic[i[0]].append(pos)     # keys = chromosomes
            else:
                dic[i[0]] = [pos]

    # Sorting tuples for each key
    for k in dic.keys():
        dic[k] = list(set(tuple(x) for x in dic[k]))
        dic[k].sort(key=lambda x: x[0])

    return dic

ref_dic = sorting_dictionary(reference_file)
print("Reference dictionary ready !")

int_dic = sorting_dictionary(interest_file)

# Testing overlap in interest dic
if int_dic.values()[1][0] != int_dic.values()[1][1]:
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
                current_ID = str(current_ID) + "," + str(new_ID)
            else:
                new_int_dic[k].append((current_start,current_end, current_ID))
                current_start = new_start
                current_end = new_end
                current_ID = new_ID

            if i == len(int_dic[k]):
                new_int_dic[k].append((current_end,current_start, current_ID))

    nb_elem = sum(len(int_dic[x]) for x in int_dic.keys())
    new_nb_elem = sum(len(new_int_dic[x]) for x in new_int_dic.keys())
    if nb_elem != new_nb_elem:
        print("There is overlap in interest file !")

    int_dic = new_int_dic

print("Reference dictionary ready !")
print("Running overlap ... ")

# Overlap interest to reference
dic_output = {}
for chr in ref_dic.keys():
    first_i = 0
    for pos in ref_dic[chr]:
        start = pos[0]
        end = pos[1]
        if chr in int_dic.keys():
            i = first_i
            ref_pos = str(chr) + "\t" + str(start) + "\t" + str(end)

            # Initialization of first possible overlapping interest position
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
                dic_output[ref_pos] = ["NA"]


output = open(output_file, 'w')
if os.stat(output_file).st_size == 0:
    output.write("chr\tstart\tend\t"+interest_name+"\n")

for ref_pos, int_pos in dic_output.items():
    output.write(ref_pos + "\t")
    for i in int_pos:
        output.write(str(i) + ",")
    output.write("\n")

output.close()
