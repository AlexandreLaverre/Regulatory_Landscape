#!/usr/bin/env python3
# coding=utf-8

import os
import re
import argparse
parser = argparse.ArgumentParser()
#parser.add_argument("specie", help="reference specie")
parser.add_argument("reference_file", help="list of genomic coordinates")
parser.add_argument("interest_file", help="list of interest genomic coordinates to be overlap")
parser.add_argument("output_file", help="name of output file")
parser.add_argument("--extend", nargs="?", default=0, const=0, type=int, help="allow for greater overlap (in base pairs) (default = 0)")
parser.add_argument("--intraoverlap", action="store_true", help="check overlap in interest file")
parser.add_argument("--count_overlap", action="store_true", help="count overlapping base pairs")
parser.add_argument("--count_window", action="store_true", help="count interest base pairs in extended windows")
parser.add_argument("--interest_ID", action="store_true", help="overlap interest output with ID format")
parser.add_argument("--reference_ID", action="store_true", help="reference output with ID format")
parser.add_argument("-v", "--verbose", action="store_true", help="increase output verbosity")
args = parser.parse_args()

#path = "/home/laverre/Documents/Regulatory_Landscape/data/" + args.specie + "/"
reference_file = args.reference_file #path +
interest_file = args.interest_file #path + 
output_file = args.output_file # path +


### Create dictionary of infiles : {chr = [(start, end, ID),...]}
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
            chr = 'chr' + str(i[0].strip('chr'))    # Check if in "chr format"
            if file == reference_file:
                ID = str(i[3]) if args.reference_ID else str(i[0])+':'+str(i[1])+':'+str(i[2])
            else:
                ID = str(i[3]) if args.interest_ID else chr + ':' + str(i[1]) + ':' + str(i[2])

            pos = (int(i[1]), int(i[2]), ID)

            if chr in dic.keys():
                dic[chr].append(pos)
            else:
                dic[chr] = [pos]

    # Sorting coordinates for each chromosome
    for k in dic.keys():
        dic[k] = list(set(tuple(x) for x in dic[k]))
        dic[k].sort(key=lambda x: x[0])

    return dic


ref_dic = sorted_dictionary(reference_file)
print("Reference dictionary ready") if args.verbose else None
int_dic = sorted_dictionary(interest_file)
print("Interest dictionary ready") if args.verbose else None


# Testing intra-overlap in interest dic
def collapse_intraoverlap(dic):
    if list(dic.values())[0][0][0] != list(dic.values())[0][0][1]:
        print("Length of", name_dic, "seq. > 1 pb ") if args.verbose else None
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
            print("There is overlap in", name_dic, "file !") if args.verbose else None
            print("Before collapse:", nb_elem, "elements; After collapse:", new_nb_elem) if args.verbose else None
        else:
            print("There is no overlap in", name_dic, "file !") if args.verbose else None

        return new_dic

    else:
        print("Length of", name_dic, "seq. = 1 pb") if args.verbose else None
        return dic


if args.intraoverlap:
    name_dic = 'interest'
    int_dic = collapse_intraoverlap(int_dic)
    print("Check intraoverlap in interest dictionary done !") if args.verbose else None

print("Running overlap... ") if args.verbose else None

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
            while i < len(int_dic[chr]) and int_dic[chr][i][1] < (start - args.extend):
                i += 1
            first_i = i

            # Adding all overlapping interest position to reference position
            while i < len(int_dic[chr]) and int_dic[chr][i][0] <= (end + args.extend):
                if ref_pos in dic_output.keys():
                    # if not any(int_dic[chr][i][2] in overlap for overlap in dic_output[ref_pos]):
                    dic_output[ref_pos].append(int_dic[chr][i])
                else:
                    dic_output[ref_pos] = [int_dic[chr][i]]
                i += 1

            # Adding reference position without overlap
            if ref_pos not in dic_output.keys():
                dic_output[ref_pos] = [('NA', 'NA', 'NA', 'NA')]

        else:
            dic_output[ref_pos] = [('NA', 'NA', 'NA', 'NA')]

        dic_output[ref_pos] = set(dic_output[ref_pos])


if args.count_overlap:
    print("Counting base pairs...")
    count_bp = {}
    length_pos = {}
    for ref_pos in dic_output.keys():
        frag_pos = ref_pos.split('\t')
        length_frag = int(frag_pos[2]) - int(frag_pos[1])
        length_pos[ref_pos] = length_frag if length_frag > 0 else 1

        for overlap in dic_output[ref_pos]:
            if overlap[2] != 'NA':
                overlap_start = int(overlap[0])
                overlap_end = int(overlap[1])

                ov_start = max(int(frag_pos[1]), overlap_start)
                ov_end = min(int(frag_pos[2]), overlap_end)
                length_overlap = ov_end - ov_start

                if ref_pos in count_bp.keys():
                    count_bp[ref_pos] = int(count_bp[ref_pos]) + length_overlap
                else:
                    count_bp[ref_pos] = length_overlap

            else:
                count_bp[ref_pos] = 0


if args.count_window:
    print("Counting base pairs in window...")
    count_nb_window = {}
    for ref_pos in dic_output.keys():
        frag_pos = ref_pos.split('\t')
        low_border = int(frag_pos[1]) - args.extend
        high_border = int(frag_pos[2]) + args.extend

        for overlap in dic_output[ref_pos]:
            if overlap[2] != 'NA':
                min_window = max(low_border, int(overlap[0]))
                max_window = min(high_border, int(overlap[1]))

                if ref_pos in count_nb_window.keys():
                    count_nb_window[ref_pos] = int(count_nb_window[ref_pos]) + (max_window-min_window)
                else:
                    count_nb_window[ref_pos] = max_window-min_window

            else:
                count_nb_window[ref_pos] = 0


print("Writting output... ") if args.verbose else None
output = open(output_file, 'w')
if os.stat(output_file).st_size == 0:
    if args.count_overlap:
        output.write("ID\tchr\tstart\tend\toverlap_ID\tlength_frag\tnb_bp_overlap\t%overlap\n") #\tgene_ID
    elif args.count_window:
        output.write("ID\tchr\tstart\tend\tnb_bp_window\n")  # \tgene_ID

    else:
        output.write("ID\tchr\tstart\tend\toverlap_ID\n")

for ref_pos, int_pos in dic_output.items():
    frag = ref_pos.split('\t')
    output.write(frag[3] + '\t' + frag[0] + "\t" + frag[1] + "\t" + frag[2] + "\t")  # chr + start + end (ref 1)
    count = 0
    gene = []
    for i in int_pos:
        count += 1
        if count == len(int_pos):
            if args.count_overlap:
                if str(i[0]) == 'NA':
                    output.write('NA' + "\t" + str(length_pos[ref_pos]) + "\t" + str(count_bp[ref_pos]) + "\t"
                                 + str((count_bp[ref_pos] / length_pos[ref_pos]) * 100) + "\n")
                else:
                    output.write(str(frag[0]) + ':' + str(i[0]) + ':' + str(i[1]) + "\t" + str(length_pos[ref_pos])
                                 + "\t" + str(count_bp[ref_pos]) + '\t' + str((count_bp[ref_pos] / length_pos[ref_pos]) * 100) + "\n")

            elif args.count_window:
                if str(i[0]) == 'NA':
                    output.write(str(count_nb_window[ref_pos]) + "\n")
                else:
                    output.write(str(count_nb_window[ref_pos]) + "\n")

            else:
                if str(i[0]) == 'NA':
                    output.write('NA' + "\n")
                else:
                    if args.interest_ID:
                        gene.append(str(i[2]))
                        gene = set(gene)
                        output.write(",".join(gene) + "\n")

                    else:
                        output.write(str(frag[0]) + ':' + str(i[0]) + ':' + str(i[1]) + "\n")  # chr:start:end

                #output.write(str(i[2]) + "\n") # only ID
                #gene.append(str(i[3]))
                #gene = set(gene)
                #output.write(str(i[2]) + "\t" + str(','.join(str(x) for x in gene)) + "\n")  # ID + gene

        else:
            #output.write(str(i[2]) + ",")
            if not args.count_window:
                if args.interest_ID:
                    gene.append(str(i[2]))
                else:
                    output.write(str(frag[0]) + ':' + str(i[0]) + ':' + str(i[1])+',')

output.close()
print("All done !") if args.verbose else None
