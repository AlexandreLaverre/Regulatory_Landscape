#!/usr/bin/env python3
# coding=utf-8

import os
import numpy as np

path_overlap = "/home/laverre/Documents/Regulatory_Landscape/data/"
path_annot = "/home/laverre/Data/Regulatory_landscape/result/Supplementary_dataset3_annotations/"
path_HIC = "/home/laverre/Data/Regulatory_landscape/result/"


def HiC_stats(origin_sp, enh):

    ############################################ Composition in base pairs ###########################################
    def dic_pb(file):
        elem_pb = {}
        with open(path_overlap + origin_sp + "/overlap/" + enh + file) as f2:
            for i in f2.readlines()[1:]:
                i = i.strip("\n")
                i = i.split("\t")
                frag = (str(i[1]) + ':' + str(i[2]) + ':' + str(i[3]))
                elem_pb[frag] = 0 if i[4] == "NA" else int(i[6])    # len(i[4].split(",")) #

        return elem_pb

    all = "_overlap_all_exons.txt"
    all_exon = dic_pb(all)

    ############################################## Duplication score ##############################################
    frag_dupli = {}
    repeat_pb = {}
    GC_pb = {}
    input = "/" + enh + "/" + enh + "_BLAT_summary_0.8.txt"

    with open(path_annot + origin_sp + input) as f1:
        for i in f1.readlines()[1:]:
            i = i.strip("\n")
            i = i.split("\t")
            frag = str(i[0])
            frag_dupli[frag] = int(i[4])
            repeat_pb[frag] = int(i[2])
            GC_pb[frag] = int(i[3])

    print("Overlap & Score dupli done ! ")

    ############################################## Contact informations ##############################################
    datas = ["_observed", "_simulated"]
    for data in datas:
        print(origin_sp, data)

        # Matching with contacted fragments
        PIR_infos = {}
        gene_infos = {}
        contact_unbaited = {}

        if data == "_simulated":
            infile = "Supplementary_dataset4_genes_enhancers_contacts/" + origin_sp + "/" + enh + "/gene_" + enh + "_enhancers_simulated_interactions.txt_unique"
        else:
            infile = "Supplementary_dataset4_genes_enhancers_contacts/" + origin_sp + "/" + enh + "/gene_" + enh + "_enhancers_original_interactions.txt_unique"

        with open(path_HIC + infile) as f3:
            first_line = f3.readline().strip("\n")
            first_line = first_line.split("\t")
            cell_name = first_line[7:]

            for i in f3.readlines():
                i = i.strip("\n")
                i = i.split("\t")
                baits = [i[2].split(",")]
                gene = str(i[0])
                PIR = str(i[1])

                dist = str(i[4])
                contact_score = [float(x) if x != "nan" else np.nan for x in i[7:]]
                contacted_gene_in_cell = [1 if x != "nan" else 0 for x in i[7:]]  #Nb gene for each cell
                median_score = str(np.median([float(x) for x in i[7:] if x != "nan"]))

                # PIR side
                if PIR not in PIR_infos.keys():
                    PIR_infos[PIR] = [baits, [dist], contact_score, contacted_gene_in_cell, [median_score], [gene]]

                else:
                    for bait in baits:
                        if bait not in PIR_infos[PIR][0]:
                            PIR_infos[PIR][0].append(bait)
                    if gene not in PIR_infos[PIR][5]:
                        PIR_infos[PIR][5].append(gene)

                    PIR_infos[PIR][1].append(str(dist))
                    PIR_infos[PIR][2] = list(np.nansum((PIR_infos[PIR][2], contact_score), axis=0)) # useless
                    PIR_infos[PIR][3] = [sum(x) for x in zip(PIR_infos[PIR][3], contacted_gene_in_cell)]
                    PIR_infos[PIR][4].append(str(median_score))

                # Bait side : names and distances of contacted fragments
                if gene not in gene_infos.keys():
                    gene_infos[gene] = [[PIR], [dist]]
                else:
                    gene_infos[gene][0].append(PIR)
                    gene_infos[gene][1].append(dist)

                # Count contacted enhancer & unbaited contact per bait
                if gene not in contact_unbaited.keys():
                    contact_unbaited[gene] = []

                if i[6] == "unbaited":
                    contact_unbaited[gene].append(PIR)

        ############################################## OUTPUT ##############################################
        print("Writting output...")
        output = open(path_annot + "/" + origin_sp + "/contacted_" + enh + "_composition_" + origin_sp + data + ".txt_unique", 'w')
        if os.stat(path_annot + "/" + origin_sp + "/contacted_" + enh + "_composition_" + origin_sp + data + ".txt_unique").st_size == 0:
            output.write("chr\tstart\tend\tlength\t")
            output.write("bait_contacted\tgenes_contacted\tmean_genes_contacts\tmedian_score\tmedian_dist\tnb_sample\t"
                         "BLAT_match\tall_exon_bp\trepeat_bp\tGC_bp\t")

            output.write('\t'.join(cell_name) + "\n")

        for PIR in PIR_infos.keys():
            nb_sample = str(len([x for x in PIR_infos[PIR][3] if float(x) > 0]))
            gene_contact = [len(gene_infos[gene][0]) for gene in PIR_infos[PIR][5]]  # nb contact of contacted gene
            length = str(int(PIR.split(':')[2]) - int(PIR.split(':')[1]))

            output.write(PIR.split(':')[0] + '\t' + PIR.split(':')[1] + '\t' + PIR.split(':')[2] + '\t' + length + '\t')
            output.write(str(len(PIR_infos[PIR][0])) + '\t' + str(len(PIR_infos[PIR][5])) + '\t' + str(np.mean(gene_contact)) + '\t')
            output.write(str(np.median([float(x) for x in PIR_infos[PIR][4]])) + '\t')  # mid_score
            output.write(str(np.median([float(x) for x in PIR_infos[PIR][1]])) + '\t')
            output.write(str(nb_sample) + '\t' + str(frag_dupli[PIR]) + '\t' + str(all_exon[PIR]) + '\t')
            output.write(str(repeat_pb[PIR]) + '\t' + str(GC_pb[PIR]) + '\t')
            output.write(str('\t'.join(str(x) for x in PIR_infos[PIR][3])) + '\n') # nb contacted gene per cell

        output.close()


origin_sp = "mouse"
enhancers = ["CAGE", "ENCODE"]
if origin_sp == "human":
    enhancers.extend(["RoadMap", "GRO_seq"])

for enh in enhancers:
    print("Running", enh)
    HiC_stats(origin_sp, enh)

print("All done ! ")
