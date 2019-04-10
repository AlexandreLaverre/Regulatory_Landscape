#!/usr/bin/env python3
# coding=utf-8

import os
import sys
import re

path_data = "/home/laverre/Documents/Regulatory_Landscape/data/human/"

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
def collapse_intraoverlap(dic):
    if list(dic.values())[0][0][0] != list(dic.values())[0][0][1]:
        print("Length of", name_dic, "seq. > 1 pb ")
        if collapse == "T":
            new_dic = {}
            for k in dic.keys():
                current_start = dic[k][0][0]
                current_end = dic[k][0][1]
                current_ID = dic[k][0][2]

                new_dic[k] = []
                if len(dic[k]) == 1:
                    new_dic[k].append((current_start, current_end, current_ID))

                for i in range(1, len(dic[k])):
                    new_start = dic[k][i][0]
                    new_end = dic[k][i][1]
                    new_ID = dic[k][i][2]
                    if current_end >= new_start >= current_start:
                        if new_end > current_end:
                            current_end = new_end
                        current_ID = str(current_ID) + "_" + str(new_ID)

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
                print("Before collapse:", nb_elem, "elements; After collapse:", new_nb_elem)
            else:
                print("There is no overlap in", name_dic, "file !")

            return new_dic

    else:
        print("Length of", name_dic, "seq. = 1 pb")


name_dic = 'reference'
ref_dic = collapse_intraoverlap(ref_dic)
name_dic = 'interest'
int_dic = collapse_intraoverlap(int_dic)


print("Interest dictionary ready")
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
                    dic_output[ref_pos].append(int_dic[chr][i])
                else:
                    dic_output[ref_pos] = [int_dic[chr][i]]
                i += 1

            # Adding reference position without overlap
            if ref_pos not in dic_output.keys():
                dic_output[ref_pos] = [('NA', 'NA', 'NA')]

        else:
            dic_output[ref_pos] = [('NA', 'NA', 'NA')]

print("Keep non overlapping regions...")

dic_no_ov = {}
for ref_pos in dic_output.keys():
    frag_pos = ref_pos.split('\t')
    chr = frag_pos[0]
    current_start = int(frag_pos[1])
    current_end = int(frag_pos[2])
    dic_no_ov[ref_pos] = []

    if dic_output[ref_pos][0][0] != "NA":
        x = 0

        while x < len(dic_output[ref_pos])-1:
            start_ov = int(dic_output[ref_pos][x][0])
            end_ov = int(dic_output[ref_pos][x][1])
            next_start = int(dic_output[ref_pos][x+1][0])

            if start_ov > current_start:
                dic_no_ov[ref_pos].append((chr + ':' + str(current_start) + ':' + str(start_ov)))
            else:
                dic_no_ov[ref_pos].append((chr + ':' + str(end_ov) + ':' + str(next_start)))

            x += 1
            current_start = end_ov

        start_ov = int(dic_output[ref_pos][x][0])
        end_ov = int(dic_output[ref_pos][x][1])

        if start_ov > current_start and end_ov < current_end:
            dic_no_ov[ref_pos].append((chr + ':' + str(current_start) + ':' + str(start_ov)))
            dic_no_ov[ref_pos].append((chr + ':' + str(end_ov) + ':' + str(current_end)))

        else:
            if start_ov < current_start:
                current_start = end_ov

            if end_ov > current_end:
                current_end = start_ov

            if current_end < current_start:
                dic_no_ov[ref_pos] = ["NULL"]
            else:
                dic_no_ov[ref_pos].append((chr + ':' + str(current_start) + ':' + str(current_end)))

    else:
        dic_no_ov[ref_pos].append((chr + ':' + str(current_start) + ':' + str(current_end)))

print("Writting output... ")

output = open(output_file, 'w')
if os.stat(output_file).st_size == 0:
    output.write("chr\tstart\tend\tno_exonic\n")

for ref_pos, int_pos in dic_no_ov.items():
    output.write(ref_pos + "\t")
    count = 0

    for i in int_pos:
        count += 1
        if count == len(int_pos):
            output.write(str(i) + "\n")
            #output.write(str(i[2])+":"+str(i[0])+"-"+str(i[1])+"\n")
        else:
            output.write(str(i) + ",")
            #output.write(str(i[2])+":"+str(i[0])+"-"+str(i[1]) + ",")

output.close()
print("All done !")
