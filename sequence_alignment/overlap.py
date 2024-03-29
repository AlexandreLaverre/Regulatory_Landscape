#!/usr/bin/env python3.7
# coding=utf-8

import os
import sys
import re


reference_file = sys.argv[1]
interest_file = sys.argv[2]
output_file = sys.argv[3]

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
            if file == reference_file:
                ID = str(i[0])+':'+str(i[1])+':'+str(i[2])+':'+str(i[5])
            else:
                ID = str(i[0])+':'+str(i[1])+':'+str(i[2]) # chr:start:end

            pos = (int(i[1]), int(i[2]), ID)  # value = position (start + end) + ID

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
print("Reference dictionary ready")

int_dic = sorted_dictionary(interest_file)


# Testing overlap in interest dic
def collapse_intraoverlap(dic):
    if list(dic.values())[0][0][0] != list(dic.values())[0][0][1]:
        print("Length of", name_dic, "seq. > 1 pb ")
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
            print("Before collapse:", nb_elem, "elements; After collapse:", new_nb_elem)
        else:
            print("There is no overlap in", name_dic, "file !")

        return new_dic

    else:
        print("Length of", name_dic, "seq. = 1 pb")
        return dic


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
        ref_pos = str(chr) + '\t' + str(start) + "\t" + str(end) + '\t' + str(pos[2])  # only ID

        if chr in int_dic.keys():

            # Initialization of first possible overlapping interest position
            i = first_i
            while i < len(int_dic[chr]) and int_dic[chr][i][1] < start:
                i += 1
            first_i = i

            # Adding all overlapping interest position to reference position
            while i < len(int_dic[chr]) and int_dic[chr][i][0] <= end:
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
    output.write("ID\tchr\tstart\tend\toverlap_ID\n")

for ref_pos, int_pos in dic_output.items():
    frag = ref_pos.split('\t')
    output.write(frag[3] + '\t' + frag[0] + "\t" + frag[1] + "\t" + frag[2] + "\t")  # chr + start + end (ref 1)
    count = 0

    for i in int_pos:
        count += 1
        if count == len(int_pos):
            if str(i[0]) == 'NA':
                output.write('NA' + "\n")
            else:
                output.write(str(frag[0]) + ':' + str(i[0]) + ':' + str(i[1]) + "\n")  # chr:start:end

        else:
            output.write(str(frag[0]) + ':' + str(i[0]) + ':' + str(i[1])+',')

output.close()
print("All done !")
