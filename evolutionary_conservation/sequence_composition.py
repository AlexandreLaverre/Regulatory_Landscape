#!/usr/bin/env python3
# coding=utf-8

import os
import numpy as np

path_overlap = "/home/laverre/Documents/Regulatory_Landscape/data/"
path_dupli = "/home/laverre/Data/Regulatory_landscape/result/Supplementary_dataset3_annotations/"
path_HIC = "/home/laverre/Data/Regulatory_landscape/result/"
path_output = "/home/laverre/Documents/Regulatory_Landscape/result/"
origin_sp = "human"


# Composition in base pairs
def dic_pb(file):
    elem_pb = {}
    with open(path_overlap + origin_sp + "/overlap/" + origin_sp + "_all_fragments" + file) as f2:
        for i in f2.readlines()[1:]:
            i = i.strip("\n")
            i = i.split("\t")
            frag = (str(i[1]) + ':' + str(i[2]) + ':' + str(i[3]))
            elem_pb[frag] = 0 if i[4] == "NA" else int(i[6])

    return elem_pb


all = "_overlap_all_exons.txt"
all_exon = dic_pb(all)
all_250 = "_overlap_all_exons250.txt"
all_exon250 = dic_pb(all_250)
coding = "_overlap_coding_exons.txt"
coding_exon = dic_pb(coding)
nocoding = "_overlap_nocoding_exons.txt"
nocoding_exon = dic_pb(nocoding)
repeat = "_overlap_repeatmasker.txt"
repeat_pb = dic_pb(repeat)
phastcons = "_overlap_phastconselements_noexonic250.txt"
phastcons_pb = dic_pb(phastcons)
CAGE = "_overlap_CAGE.txt"
CAGE_count = dic_pb(CAGE)
ENCODE = "_overlap_ENCODE.txt"
ENCODE_count = dic_pb(ENCODE)
if origin_sp == "human":
    RoadMap = "_overlap_RoadMap.txt"
    RoadMap_count = dic_pb(RoadMap)
    GRO_seq = "_overlap_GRO_seq.txt"
    GRO_seq_count = dic_pb(GRO_seq)


# Composition in count number
def dic_count(file):
    elem_pb = {}
    with open(path_overlap + origin_sp + "/overlap/" + origin_sp + "_all_fragments" + file) as f2:
        for i in f2.readlines()[1:]:
            i = i.strip("\n")
            i = i.split("\t")
            frag = (str(i[1]) + ':' + str(i[2]) + ':' + str(i[3]))
            if file == "_overlap_genes_1Kb.txt":
                elem_pb[frag] = [] if i[4] == "NA" else i[4].split(",")
            else:
                elem_pb[frag] = 0 if i[4] == "NA" else len(i[4].split(","))

    return elem_pb


TSS = "_overlap_genes.txt"
TSS_count = dic_count(TSS)
bait_TSS1Kb = "_overlap_TSS_1Kb.txt"
TSS_count_bait = dic_count(bait_TSS1Kb)
bait_gene1Kb = "_overlap_genes_1Kb.txt"
gene_count_bait = dic_count(bait_gene1Kb)

# Duplication score
frag_dupli = {}
input = "/all_fragments_BLAT_summary_0.8.txt"

with open(path_dupli + origin_sp + input) as f1:
    for i in f1.readlines()[1:]:
        i = i.strip("\n")
        i = i.split("\t")
        frag = i[0].split(':')
        frag = str(frag[0]) + ':' + str(frag[1]) + ':' + str(frag[2])
        frag_dupli[frag] = int(i[4]) - 1

print("Overlap & Score dupli done ! ")


def HiC_stats(origin_sp, data):
    print(origin_sp, data)

    # Matching with contacted fragments
    PIR_infos = {}
    bait_infos = {}
    contact_unbaited = {}
    CAGE_contact = {}
    RoadMap_contact = {}
    ENCODE_contact = {}
    GRO_seq_contact = {}

    if data == "_simulated":
        infile = "Supplementary_dataset2_simulated_interactions/" + origin_sp + "/" + "simulated_all_interactions.txt"
    else:
        infile = "Supplementary_dataset1_original_interactions/" + origin_sp + "/" + "all_interactions.txt"

    with open(path_HIC + infile) as f3:
        first_line = f3.readline().strip("\n")
        first_line = first_line.split("\t")
        cell_name = first_line[8:]

        for i in f3.readlines():
            i = i.strip("\n")
            i = i.split("\t")
            bait = (str(i[0]) + ":" + str(i[1]) + ":" + str(i[2]))
            PIR = (str(i[3]) + ":" + str(i[4]) + ":" + str(i[5]))

            if i[0] == i[3]:
                dist = float(i[7])
                if 25000 <= dist <= 10000000:
                    contact = [float(x) if x != "NA" else np.nan for x in i[8:]]
                    nb_bait_in_cell = [1 if x != "NA" else 0 for x in i[8:]]  #  Nb bait for each cell
                    score = str(np.median([float(x) for x in i[8:] if x != "NA"]))

                    # PIR side
                    if PIR not in PIR_infos.keys():
                        PIR_infos[PIR] = [[bait], [str(dist)], contact, nb_bait_in_cell, i[6], [score]]
                    else:
                        PIR_infos[PIR][0].append(bait)
                        PIR_infos[PIR][1].append(str(dist))
                        PIR_infos[PIR][2] = list(np.nansum((PIR_infos[PIR][2], contact), axis=0))
                        PIR_infos[PIR][3] = [sum(x) for x in zip(PIR_infos[PIR][3], nb_bait_in_cell)]
                        PIR_infos[PIR][5].append(str(score))

                    # Bait side : names and distances of contacted fragments
                    if bait not in bait_infos.keys():
                        bait_infos[bait] = [[PIR], [str(dist)]]
                    else:
                        bait_infos[bait][0].append(PIR)
                        bait_infos[bait][1].append(str(dist))

                    # Count contacted enhancer & unbaited contact per bait
                    if bait not in contact_unbaited.keys():
                        contact_unbaited[bait] = []
                        CAGE_contact[bait] = 0
                        RoadMap_contact[bait] = 0
                        ENCODE_contact[bait] = 0
                        GRO_seq_contact[bait] = 0

                    if i[6] == "unbaited":
                        contact_unbaited[bait].append(PIR)

                    CAGE_contact[bait] += CAGE_count[PIR]
                    ENCODE_contact[bait] += ENCODE_count[PIR]

                    if origin_sp == "human":
                        RoadMap_contact[bait] += RoadMap_count[PIR]
                        GRO_seq_contact[bait] += GRO_seq_count[PIR]

    print("Writting output...")
    output = open("../../result/conservation/contacted_sequence_composition_" + origin_sp + data + ".txt_test", 'w')
    if os.stat("../../result/conservation/contacted_sequence_composition_" + origin_sp + data + ".txt_test").st_size == 0:
        output.write("chr\tstart\tend\tlength\tCAGE_bp\tENCODE_bp\t")
        if origin_sp == "human":
            output.write("RoadMap_bp\tGRO_seq_bp\t")
        output.write("bait_contacted\tmean_baits_contacts\tmedian_score\tmidist\tbaited\tnb_cell\tduplication\tall_exon_pb\t"
                     "all_exon250\tcoding_exon_pb\tnocoding_exon_pb\trepeat_pb\tphastcons_noexonic250\tTSS_count\t")
        output.write('\t'.join(cell_name) + "\n")

    for PIR in PIR_infos.keys():
        nb_sample = str(len([x for x in PIR_infos[PIR][2] if float(x) > 0]))
        bait_contact = [len(bait_infos[bait][0]) for bait in PIR_infos[PIR][0]]  # nb contact of contacted baits
        length = str(int(PIR.split(':')[2]) - int(PIR.split(':')[1]))

        output.write(PIR.split(':')[0] + '\t' + PIR.split(':')[1] + '\t' + PIR.split(':')[2] + '\t' + length + '\t')
        output.write(str(CAGE_count[PIR]) + '\t' + str(ENCODE_count[PIR]) + '\t')
        if origin_sp == "human":
            output.write(str(RoadMap_count[PIR]) + '\t' + str(GRO_seq_count[PIR]) + '\t')

        output.write(str(len(PIR_infos[PIR][0])) + '\t' + str(np.mean(bait_contact)) + '\t')
        output.write(str(np.median([float(x) for x in PIR_infos[PIR][5]])) + '\t')  # mid_score
        output.write(str(np.median([float(x) for x in PIR_infos[PIR][1]])) + '\t' + str(PIR_infos[PIR][4]) + '\t')
        output.write(str(nb_sample) + '\t' + str(frag_dupli[PIR]) + '\t' + str(all_exon[PIR]) + '\t')
        output.write(str(all_exon250[PIR]) + '\t' + str(coding_exon[PIR]) + '\t' + str(nocoding_exon[PIR]) + '\t')
        output.write(str(repeat_pb[PIR]) + '\t' + str(phastcons_pb[PIR]) + '\t' + str(TSS_count[PIR]) + '\t')
        output.write(str('\t'.join(str(x) for x in PIR_infos[PIR][3])) + '\n')

    ## Bait side
    output_bait = open("../../result/conservation/bait_composition_" + origin_sp + data + ".txt_test", 'w')
    if os.stat("../../result/conservation/bait_composition_" + origin_sp + data + ".txt_test").st_size == 0:
        output_bait.write("chr\tstart\tend\tCAGE_contacted\tENCODE_contacted\t")
        if origin_sp == "human":
            output_bait.write("RoadMap_contacted\tGRO_seq_contacted\t")
        output_bait.write("PIR_contacted\tunbaited_PIR_contacted\tPIR_complexity\tmidist\tduplication\tall_exon_pb"
                          "\tcoding_exon_pb\tnocoding_exon_pb\trepeat_pb\tphastcons_noexonic250\tTSS_count1Kb"
                          "\tgenes_count1Kb\tgenes\n")

    for Bait in bait_infos.keys():
        PIR_contact = [len(PIR_infos[PIR][0]) for PIR in bait_infos[Bait][0]]  # nb contact of contacted PIR

        output_bait.write(Bait.split(':')[0] + '\t' + Bait.split(':')[1] + '\t' + Bait.split(':')[2]+ '\t')
        output_bait.write(str(CAGE_contact[Bait]) + '\t' + str(ENCODE_contact[Bait]) + '\t')
        if origin_sp == "human":
            output_bait.write(str(RoadMap_contact[Bait]) + '\t' + str(GRO_seq_count[Bait]) + '\t')

        output_bait.write(str(len(bait_infos[Bait][0])) + '\t' + str(len(contact_unbaited[Bait])) + '\t' +
                          str(np.mean(PIR_contact)) + '\t' + str(np.median([float(x) for x in bait_infos[Bait][1]]))
                          + '\t' + str(frag_dupli[Bait]) + '\t' + str(all_exon[Bait]) + '\t' + str(all_exon250[Bait])
                          + '\t' + str(coding_exon[Bait]) + '\t' + str(nocoding_exon[Bait])+ '\t' + str(repeat_pb[Bait])
                          + '\t' + str(phastcons_pb[Bait]) + '\t' + str(TSS_count_bait[Bait])
                          + '\t' + str(len(gene_count_bait[Bait])) + '\t' + str(",".join(gene_count_bait[Bait])) + '\n')

    output.close()


datas = ["_observed", "_simulated"]

for data in datas:
    HiC_stats(origin_sp, data)

print("All done ! ")
