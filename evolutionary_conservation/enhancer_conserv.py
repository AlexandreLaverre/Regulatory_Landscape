#!/usr/bin/env python3
# coding=utf-8

import os
from collections import defaultdict
from itertools import chain

path = "/home/laverre/Data/Regulatory_landscape/result"
path_evol = path + "/Supplementary_dataset6_regulatory_landscape_evolution/"
path_annot = path + "/Supplementary_dataset3_annotations/"


############################################# Enhancers alignments ###################################################
def enh_score(enh_name, target_sp):
    align = {}
    with open(path_evol + ref_sp + "/enhancers_conservation/" + enh_name +
              "/AlignmentStatistics_Excluding_all_Exons_" + ref_sp + "2" + target_sp + "_" + enh_name + ".txt") as f1:
        for i in f1.readlines()[1:]:
            i = i.strip("\n")
            i = i.split("\t")

            enh_origin = i[0].strip(":+")
            try:
                align_score = int(i[6]) / int(i[9])  # FilteredUngappedLength / FilteredAlignmentLength
            except ZeroDivisionError:
                align_score = 0

            align[enh_origin] = str(align_score)

    return align


def enh_conserv(enh_name):
    # Get alignment score for all species
    align_cow = enh_score(enh_name, "cow")
    align_opossum = enh_score(enh_name, "opossum")
    align_elephant = enh_score(enh_name, "elephant")
    align_rabbit = enh_score(enh_name, "rabbit")
    align_rat = enh_score(enh_name, "rat")
    align_macaque = enh_score(enh_name, "macaque")
    align_dog = enh_score(enh_name, "dog")
    align_chicken = enh_score(enh_name, "chicken")
    align_other = enh_score(enh_name, "mouse") if ref_sp == "human" else enh_score(enh_name, "human")

    all = [align_other, align_macaque, align_cow, align_elephant, align_rabbit, align_rat, align_dog, align_opossum, align_chicken]

    # Get all enhancers keys
    all_enh = {}
    with open(path_annot + ref_sp + "/" + enh_name + "_BLAT_summary_0.8.txt") as f1:
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
    for enh, score in chain(align_other.items(), align_macaque.items(), align_cow.items(), align_elephant.items(),
                            align_rabbit.items(), align_rat.items(), align_dog.items(), align_opossum.items(),
                            align_chicken.items()):

        align_all_sp[enh].append(score)

    output_file = path_evol + ref_sp + "/enhancers_conservation/" + enh_name + "/Alignments_stats_all_species.txt"
    output = open(output_file, 'w')
    if os.stat(output_file).st_size == 0:
        other = "human" if ref_sp == "mouse" else "mouse"
        output.write("enh\t"+other+"\tmacaque\tcow\telephant\trabbit\trat\tdog\topossum\tchicken\n")

    for enh, scores in align_all_sp.items():
        output.write(enh + '\t' + "\t".join(scores) + '\n')

    output.close()


references_sp = ["mouse", "human"]
enh_datas = ["CAGE"]

for ref_sp in references_sp:
    if ref_sp == "human":
        enh_datas.extend(["ENCODE", "GRO_seq", "RoadMap"])

    for enh_data in enh_datas:
        enh_conserv(enh_data)
        print(enh_data, "in", ref_sp, "done !")

