#!/usr/bin/env python3
# coding=utf-8

import os
import sys
import re

path_data = "/home/laverre/Documents/Regulatory_Landscape/data/mouse/"

reference_file = path_data + sys.argv[1]
interest_file = path_data + sys.argv[2]
output_file = path_data + sys.argv[3]
intraoverlap = sys.argv[4]
counting = "F"

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
            if file == reference_file:
                ID = str(i[0])+':'+str(i[1])+':'+str(i[2])+':'+str(i[5])
                #gene = "pb"
            else:
                ID = chr + ':' + str(i[1]) + ':' + str(i[2])  # chr:start:end
                #ID = str(i[3])  # only ID
                #gene = str(i[4])

            pos = (int(i[1]), int(i[2]), ID) #, gene)  # value = position (start + end) + ID

            if chr in dic.keys():
                dic[chr].append(pos)     # keys = chromosome
            else:
                dic[chr] = [pos]

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
            print("Before collapse:", nb_elem, "elements; After collapse:", new_nb_elem)
        else:
            print("There is no overlap in", name_dic, "file !")

        return new_dic

    else:
        print("Length of", name_dic, "seq. = 1 pb")
        return dic


if intraoverlap == "T":
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
            while i < len(int_dic[chr]) and int_dic[chr][i][1] < start:  # -1000
                i += 1
            first_i = i

            # Adding all overlapping interest position to reference position
            while i < len(int_dic[chr]) and int_dic[chr][i][0] <= end:  # +1000
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
                dic_output[ref_pos] = [('NA', 'NA', 'NA', 'NA')]

        else:
            dic_output[ref_pos] = [('NA', 'NA', 'NA', 'NA')]

if counting == "T":
    print("Counting base pair...")
    count_bp = {}
    length_pos = {}
    for ref_pos in dic_output.keys():
        for overlap in dic_output[ref_pos]:
            frag_pos = ref_pos.split('\t')
            length_frag = int(frag_pos[2]) - int(frag_pos[1])
            length_pos[ref_pos] = length_frag

            if length_pos[ref_pos] == 0:
                length_pos[ref_pos] = 1

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
    output.write("ID\tchr\tstart\tend\toverlap_ID\n") #\tlength_frag\tnb_bp_overlap\t%TSS_ID\tgene_ID

for ref_pos, int_pos in dic_output.items():
    frag = ref_pos.split('\t')
    output.write(frag[3] + '\t' + frag[0] + "\t" + frag[1] + "\t" + frag[2] + "\t")  # chr + start + end (ref 1)
    count = 0
    gene = []
    for i in int_pos:
        count += 1
        if count == len(int_pos):
            if str(i[0]) == 'NA':
                output.write('NA' + "\n")
            else:
                output.write(str(frag[0]) + ':' + str(i[0]) + ':' + str(i[1]) + "\n") # chr:start:end
                #output.write(str(i[2]) + "\n") # only ID
                #gene.append(str(i[3]))
                #gene = set(gene)
                #output.write(str(i[2]) + "\t" + str(','.join(str(x) for x in gene)) + "\n")  # ID + gene

                #output.write(str(frag[0]) + ':' + str(i[0]) + ':' + str(i[1]) + "\t" + str(length_pos[ref_pos]) + "\t" +
                #             str(count_bp[ref_pos]) + '\t' + str((count_bp[ref_pos] / length_pos[ref_pos]) * 100) + "\n")

        else:
            #gene.append(str(i[3]))
            #output.write(str(i[2]) + ",")
            output.write(str(frag[0]) + ':' + str(i[0]) + ':' + str(i[1])+',')

output.close()
print("All done !")
