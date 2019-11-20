#!/usr/bin/env python3
# coding=utf-8

import os
import sys
import re

path_data = "/home/laverre/Documents/Regulatory_Landscape/data/mouse/"

interest_file = path_data + sys.argv[1]
output_file = path_data + sys.argv[2]

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
            chr = 'chr' + str(i[0].strip('chr'))
            ID = chr + ':' + str(i[1]) + ':' + str(i[2])  # chr:start:end

            pos = (int(i[1]), int(i[2]), ID) # value = position (start + end) + ID

            if chr in dic.keys():
                dic[chr].append(pos)     # keys = chromosome
            else:
                dic[chr] = [pos]

    # Sorting tuples for each key
    for k in dic.keys():
        dic[k] = list(set(tuple(x) for x in dic[k]))
        dic[k].sort(key=lambda x: x[0])

    return dic


# Testing overlap in interest dic
def collapse_intraoverlap(dic):
    if list(dic.values())[0][0][0] != list(dic.values())[0][0][1]:
        print("Length of", name_dic, "seq. > 1 pb ")
        new_dic = {}
        for k in dic.keys():
            current_start = dic[k][0][0]  # - 250
            current_end = dic[k][0][1]  # + 250
            current_ID = dic[k][0][2]

            new_dic[k] = []
            if len(dic[k]) == 1:
                new_dic[k].append((current_start, current_end, current_ID))

            for i in range(1, len(dic[k])):
                new_start = dic[k][i][0]  # - 250
                new_end = dic[k][i][1]  # + 250
                new_ID = dic[k][i][2]
                if current_end >= new_start >= current_start:
                    current_ID = str(current_ID) + "_" + str(new_ID)

                    if new_end > current_end:
                        current_end = new_end

                else:
                    new_dic[k].append((current_start, current_end, current_ID))
                    current_start = new_start
                    current_end = new_end
                    current_ID = new_ID

                if i == len(dic[k])-1:
                    new_dic[k].append((current_start, current_end, current_ID))

        nb_elem = sum(len(dic[x]) for x in dic.keys())
        new_nb_elem = sum(len(new_dic[x]) for x in new_dic.keys())
        if nb_elem != new_nb_elem:
            print("There is overlap in", name_dic, "file !")
            print("Before collapse:", nb_elem, "elements")
            print("After collapse:", new_nb_elem)
        else:
            print("There is no overlap in", name_dic, "file !")

        return new_dic

    else:
        print("Length of", name_dic, "seq. = 1 pb")
        return dic


int_dic = sorted_dictionary(interest_file)
name_dic = 'interest'
int_dic = collapse_intraoverlap(int_dic)

print("Writting output... ")

output = open(output_file, 'w')

for chr in int_dic.keys():
    for frag in int_dic[chr]:
        ID = str(chr) + ":" + str(frag[0]) + ":" + str(frag[1])
        output.write(str(chr) + "\t" + str(frag[0]) + "\t" + str(frag[1]) + "\t" + ID + "\t" + "+" + "\t" + "+" + '\n')

output.close()


