#!/usr/bin/env python3
# coding=utf-8

import os
from collections import defaultdict
from itertools import chain

path = "/home/laverre/Data/Regulatory_landscape/result"
path_evol = path + "/Supplementary_dataset6_regulatory_landscape_evolution/"
path_annot = path + "/Supplementary_dataset3_annotations/"


############################################# Enhancers alignments ###################################################
def enh_score(enh_name, target_sp, score):
    align = {}
    with open(path_evol + ref_sp + "/sequence_conservation/" + enh_name +
              "/AlignmentStatistics_Excluding_all_Exons_" + ref_sp + "2" + target_sp + "_" + enh_name + ".txt") as f1:
        for i in f1.readlines()[1:]:
            i = i.strip("\n")
            i = i.split("\t")

            enh_origin = i[0].strip(":+")
            try:
                if score == "total_ungapped":
                    align_score = int(i[4]) / int(i[8])  # TotalUngappedLength / TotalAlignmentLength
                elif score == "total_identical":
                    align_score = int(i[5]) / int(i[8])  # TotalIdenticalLength / TotalAlignmentLength
                elif score == "filtered_ungapped":
                    align_score = int(i[6]) / int(i[9])  # FilteredUngappedLength / FilteredAlignmentLength
                elif score == "filtered_identical":
                    align_score = int(i[7]) / int(i[9])  # FilteredIdenticalLength / FilteredAlignmentLength

            except ZeroDivisionError:
                align_score = 0

            align[enh_origin] = str(align_score)

    return align


def enh_conserv(enh_name, score):
    # Get alignment score for all species
    align_cow = enh_score(enh_name, "cow", score)
    align_opossum = enh_score(enh_name, "opossum", score)
    align_elephant = enh_score(enh_name, "elephant", score)
    align_rabbit = enh_score(enh_name, "rabbit", score)
    align_rat = enh_score(enh_name, "rat", score)
    align_macaque = enh_score(enh_name, "macaque", score)
    align_dog = enh_score(enh_name, "dog", score)
    align_chicken = enh_score(enh_name, "chicken", score)
    align_other = enh_score(enh_name, "mouse", score) if ref_sp == "human" else enh_score(enh_name, "human", score)

    all = [align_other, align_macaque, align_cow, align_elephant, align_rabbit, align_rat, align_dog, align_opossum, align_chicken]

    # Get all sequence keys
    all_enh = {}
    with open(path_annot + ref_sp + "/" + enh_name + "/" + enh_name + "_BLAT_summary_0.8.txt") as f1:
        for i in f1.readlines()[1:]:
            i = i.strip("\n")
            i = i.split("\t")
            all_enh[i[0]] = "key"

    for dic in all:
        for enh_key in all_enh.keys():
            if enh_key not in dic:
                dic[enh_key] = str(0)

    ### Fusion of all align dictionaries
    align_all_sp = defaultdict(list)
    for enh, align_score in chain(align_other.items(), align_macaque.items(), align_cow.items(), align_elephant.items(),
                            align_rabbit.items(), align_rat.items(), align_dog.items(), align_opossum.items(),
                            align_chicken.items()):

        align_all_sp[enh].append(align_score)

    output_file = path_evol + ref_sp + "/sequence_conservation/" + enh_name + "/Alignments_stats_all_species_" + score + ".txt3"
    output = open(output_file, 'w')
    if os.stat(output_file).st_size == 0:
        other = "human" if ref_sp == "mouse" else "mouse"
        output.write("enh\t"+other+"\tmacaque\tcow\telephant\trabbit\trat\tdog\topossum\tchicken\n")

    for enh, align_scores in align_all_sp.items():
        output.write(enh + '\t' + "\t".join(align_scores) + '\n')

    output.close()


references_sp = ["mouse"] #, "human"]
enh_datas = ["ENCODE"] # "restriction_fragments", "CAGE",
scores = ["total_ungapped", "total_identical", "filtered_ungapped", "filtered_identical"]

for ref_sp in references_sp:
    #if ref_sp == "human":
    #    enh_datas.extend(["GRO_seq", "RoadMap"])

    for enh_data in enh_datas:
        for score in scores:
            enh_conserv(enh_data, score)
            print(enh_data, "in", ref_sp, "for", score, "done !")

