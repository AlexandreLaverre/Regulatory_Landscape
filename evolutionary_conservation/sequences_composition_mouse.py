#!/usr/bin/env python3
# coding=utf-8

import os
import numpy as np


def composition_seq(origin_sp, data):
    print(origin_sp, data)

    # Composition in base pairs
    def dic_pb(file):
        elem_pb = {}
        with open("../../data/"+origin_sp+"/overlap/"+origin_sp+file) as f2:
            for i in f2.readlines()[1:]:
                i = i.strip("\n")
                i = i.split("\t")
                frag = (str(i[1]) + ':' + str(i[2]) + ':' + str(i[3]))
                elem_pb[frag] = 0 if i[4] == "NA" else int(i[6])

        return elem_pb

    all = "_merged_overlap_all_exons_bp.txt"
    all_exon = dic_pb(all)
    all_250 = "_merged_overlap_all_exons250_bp.txt"
    all_exon250 = dic_pb(all_250)
    coding = "_merged_overlap_coding_exons_bp.txt"
    coding_exon = dic_pb(coding)
    nocoding = "_merged_overlap_nocoding_exons_bp.txt"
    nocoding_exon = dic_pb(nocoding)
    repeat = "_merged_overlap_repeatmasker_bp.txt"
    repeat_pb = dic_pb(repeat)
    phastcons = "_merged_overlap_phastcons_noexonic250_bp.txt"
    phastcons_pb = dic_pb(phastcons)

    # Composition in count number
    def dic_count(file):
        elem_pb = {}
        with open("../../data/"+origin_sp+"/overlap/"+origin_sp+file) as f2:
            for i in f2.readlines()[1:]:
                i = i.strip("\n")
                i = i.split("\t")
                frag = (str(i[1]) + ':' + str(i[2]) + ':' + str(i[3]))
                if file == "_bait_overlap_genes_1Kb.txt":
                    elem_pb[frag] = [] if i[4] == "NA" else i[4].split(",")
                else:
                    elem_pb[frag] = 0 if i[4] == "NA" else len(i[4].split(","))
        return elem_pb

    TSS = "_merged_overlap_genes.txt"
    TSS_count = dic_count(TSS)
    bait_TSS1Kb = "_bait_overlap_TSS_1Kb.txt"
    TSS_count_bait = dic_count(bait_TSS1Kb)
    bait_gene1Kb = "_bait_overlap_genes_1Kb.txt"
    gene_count_bait = dic_count(bait_gene1Kb)
    CAGE = "_merged_overlap_CAGE.txt"
    CAGE_count = dic_count(CAGE)

    print("Calcul overlap done ! ")

    # Duplication score
    frag_dupli = {}
    with open("../../result/BLAT_duplication/" + origin_sp + "_merged_fragments_duplication_stats.txt") as f1:
        for i in f1.readlines()[1:]:
            i = i.strip("\n")
            i = i.split("\t")
            frag = i[0].split(':')
            frag = str(frag[0]) + ':' + str(frag[1]) + ':' + str(frag[2])
            frag_dupli[frag] = int(i[3]) - 1

    print("Score dupli done ! ")

    # Matching with contacted fragments
    inter = {}
    inter_bait = {}
    link = {}
    link_bait = {}
    contact_unbaited = {}
    PIR_baited = {}
    PIR_cell = {}
    PIR_nbcell = {}
    CAGE_contact = {}

    if data == "_simul":
        infile = "/Simulations/simulations_"+origin_sp+"_10Mb_bin5kb_fragoverbin_chr"+merge+".txt"
    else:
        infile = "/all_interactions/all_interactions_chr"+merge+".txt_cell_names"

    with open("../../data/"+origin_sp+infile) as f3:

        first_line = f3.readline().strip("\n")
        first_line = first_line.split("\t")
        cell_name = first_line[13:len(first_line)]

        for i in f3.readlines():
            i = i.strip("\n")
            i = i.split("\t")
            bait = (str(i[0]) + ":" + str(i[1]) + ":" + str(i[2]))
            merged_PIR = (str(i[3]) + ":" + str(i[4]) + ":" + str(i[5]))

            if i[0] == i[3]:
                dist_obs = float(i[7])
                if 25000 <= dist_obs <= 10000000:
                    PIR_baited[merged_PIR] = i[6]

                    # Union des types cellulaires dans lesquels le contact est observé pour au mois 1 bait
                    if merged_PIR in PIR_nbcell.keys():
                        PIR_nbcell[merged_PIR] = list(set(PIR_nbcell[merged_PIR] + i[12].split(',')))
                    else:
                        PIR_nbcell[merged_PIR] = i[12].split(',')

                    # Nombre de bait du contact pour chaque cellule
                    if data == "":
                        if merged_PIR in PIR_cell.keys():
                            PIR_cell[merged_PIR] = [sum(x) for x in zip(PIR_cell[merged_PIR], list(map(int, i[13:len(i)])))]
                        else:
                            PIR_cell[merged_PIR] = list(map(int, i[13:len(i)]))
                    else:
                        PIR_cell[merged_PIR] = ["NA"]
                        cell_name = ["NA"]

                    # Bait side : names and distances of contacted fragments
                    if bait not in link_bait.keys():
                        link_bait[bait] = [merged_PIR]
                    else:
                        link_bait[bait].append(merged_PIR)
                    if bait not in inter_bait.keys():
                        inter_bait[bait] = [dist_obs]
                    else:
                        inter_bait[bait].append(dist_obs)

                    # PIR side : names and distances of target baits
                    if merged_PIR not in link.keys():
                        link[merged_PIR] = [bait]
                    else:
                        link[merged_PIR].append(bait)
                    if merged_PIR not in inter.keys():
                        inter[merged_PIR] = [dist_obs]
                    else:
                        inter[merged_PIR].append(dist_obs)

                    # Adding non-dupli fragments
                    if merged_PIR not in frag_dupli.keys():
                        frag_dupli[merged_PIR] = 0
                    if bait not in frag_dupli.keys():
                        frag_dupli[bait] = 0

                    # Count contacted enhancer & unbaited contact per bait
                    if bait not in contact_unbaited.keys():
                        contact_unbaited[bait] = []
                        CAGE_contact[bait] = 0

                    if i[6] == "unbaited":
                        contact_unbaited[bait].append(merged_PIR)

                    CAGE_contact[bait] += CAGE_count[merged_PIR]


    print("Writting output...")
    output = open("../../result/conservation/contacted_sequence_composition_"+origin_sp+data+merge+".txt", 'w')
    if os.stat("../../result/conservation/contacted_sequence_composition_"+origin_sp+data+merge+".txt").st_size == 0:
        output.write("chr\tstart\tend\tCAGE_count\tbait_contacted\tbait_complexity\tmidist_obs\tbaited\tnb_cell"
                     "\tduplication\tall_exon_pb\tall_exon250\tcoding_exon_pb\tnocoding_exon_pb\trepeat_pb"
                     "\tphastcons_noexonic250\tTSS_count\t"+'\t'.join(cell_name)+"\n")

    for PIR in inter.keys():
        bait_contact = [len(inter_bait[bait]) for bait in link[PIR]]

        output.write(PIR.split(':')[0] + '\t' + PIR.split(':')[1] + '\t' + PIR.split(':')[2]
                     + '\t' + str(CAGE_count[PIR]) + '\t' + str(len(inter[PIR]))
                     + '\t' + str(np.mean(bait_contact)) + '\t' + str(np.median(inter[PIR]))
                     + '\t' + str(PIR_baited[PIR]) + '\t' + str(len(PIR_nbcell[PIR]))
                     + '\t' + str(frag_dupli[PIR]) + '\t' + str(all_exon[PIR]) + '\t' + str(all_exon250[PIR])
                     + '\t' + str(coding_exon[PIR]) + '\t' + str(nocoding_exon[PIR]) + '\t' + str(repeat_pb[PIR])
                     + '\t' + str(phastcons_pb[PIR]) + '\t' + str(TSS_count[PIR])
                     + '\t' + str('\t'.join(str(x) for x in PIR_cell[PIR])) + '\n')

    ## Bait side
    output_bait = open("../../result/conservation/bait_composition_" + origin_sp + data + merge+".txt", 'w')
    if os.stat("../../result/conservation/bait_composition_"+origin_sp+data+merge+".txt").st_size == 0:
        output_bait.write("chr\tstart\tend\tCAGE_contacted\tPIR_contacted\tunbaited_PIR_contacted\tPIR_complexity"
                          "\tmidist_obs\tduplication\tall_exon_pb\tall_exon250\tcoding_exon_pb\tnocoding_exon_pb\trepeat_pb"
                          "\tphastcons_noexonic250\tTSS_count\tgenes_count1Kb\tgenes\n")

    for Bait in inter_bait.keys():
        PIR_contact = [len(inter[PIR]) for PIR in link_bait[Bait]]
        output_bait.write(Bait.split(':')[0] + '\t' + Bait.split(':')[1] + '\t' + Bait.split(':')[2]
                          + '\t' + str(CAGE_contact[Bait]) + '\t' + str(len(inter_bait[Bait]))
                          + '\t' + str(len(contact_unbaited[Bait])) + '\t' + str(np.mean(PIR_contact))
                          + '\t' + str(np.median(inter_bait[Bait])) + '\t' + str(frag_dupli[Bait])
                          + '\t' + str(all_exon[Bait]) + '\t' + str(all_exon250[Bait]) + '\t' + str(coding_exon[Bait])
                          + '\t' + str(nocoding_exon[Bait]) + '\t' + str(repeat_pb[Bait])
                          + '\t' + str(phastcons_pb[Bait]) + '\t' + str(TSS_count_bait[Bait])
                          + '\t' + str(len(gene_count_bait[Bait])) + '\t' + str(",".join(gene_count_bait[Bait])) + '\n')

    output.close()


origin_sp = "mouse"
datas = ["", "_simul"]
merge = "_merged"

for data in datas:
    composition_seq("mouse", data)

print("All done ! ")
