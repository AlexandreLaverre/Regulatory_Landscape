#!/usr/bin/env python3
# coding=utf-8

import os
import numpy as np

ref_sp = "human"
target_sp = "mouse" if ref_sp == "human" else "human"
seuil = 0.4

path = "/home/laverre/Data/Regulatory_landscape/result"
path_evol = path + "/Supplementary_dataset6_regulatory_landscape_evolution/" + ref_sp + "/"
path_annot = path + "/Supplementary_dataset3_annotations/" + ref_sp + "/"
path_contact = path + "/Supplementary_dataset4_genes_enhancers_contacts/" + ref_sp + "/"

########################################### Orthologous gene one2one ###############################################
ortho = {}
with open(path_evol + "gene_orthology/" + ref_sp + "2" + target_sp + "_orthologue.txt") as f3:
    for i in f3.readlines()[1:]:
        i = i.strip("\n")
        i = i.split("\t")
        if i[9] == "ortholog_one2one":
            ortho[str(i[0])] = str(i[5])


############################################# Enhancers alignments ###################################################
def enh_score(enh_name):
    align = {}
    with open(path_evol + "enhancers_conservation/" + enh_name +
              "/AlignmentStatistics_Excluding_all_Exons_" + ref_sp + "2" + target_sp + "_" + enh_name + ".txt") as f1:
        for i in f1.readlines()[1:]:
            i = i.strip("\n")
            i = i.split("\t")

            enh_origin = i[0].strip(":+")
            try:
                align_score = int(i[6]) / int(i[9])  # FilteredUngappedLength / FilteredAlignmentLength
            except ZeroDivisionError:
                align_score = 0

            if align_score > seuil:
                align[enh_origin] = align_score

    return align


############################################# Synteny conservation ###################################################
def enh_synteny(enh_name, data):
    synteny_10M = {}
    synteny_2M = {}
    with open(path_evol + "synteny_conservation/" + enh_name + "/" + ref_sp + "2" + target_sp + "_"
              + enh_name + "_" + data + "_synteny.txt") as f1:
        for i in f1.readlines()[1:]:
            i = i.strip("\n")
            i = i.split("\t")

            if float(i[12]) > seuil:
                if i[7] != "trans" and float(i[7]) < 10000000:
                    gene = i[0]
                    enh = i[2]
                    if gene not in synteny_10M.keys():
                        synteny_10M[gene] = [enh]
                    else:
                        synteny_10M[gene].append(enh)

                    if float(i[7]) < 2000000:
                        if gene not in synteny_2M.keys():
                            synteny_2M[gene] = [enh]
                        else:
                            synteny_2M[gene].append(enh)

    return synteny_10M, synteny_2M


############################################# Contact conservation ###################################################
def enh_contact(enh_name, data):
    contact = {}
    with open(path_evol + "contact_conservation/" + enh_name + "/" + ref_sp + "2" + target_sp +
              "_" + data + ".txt") as f1:
        for i in f1.readlines()[1:]:
            i = i.strip("\n")
            i = i.split("\t")

            if float(i[11]) > seuil:
                gene = i[0]
                enh = i[1]
                if gene not in contact.keys():
                    contact[gene] = [enh]
                else:
                    contact[gene].append(enh)

    return contact


############################################# Nb contact ###################################################
def summary(enh_name, data):
    align_score = enh_score(enh_name)
    enh_synt10M, enh_synt2M = enh_synteny(enh_name, data)
    enh_cont = enh_contact(enh_name, data)
    enh_total = {}
    enh_conserv = {}
    with open(path_contact + "gene_" + enh_name + "_enhancers_" + data + "_interactions.txt") as f1:
        for i in f1.readlines()[1:]:
            i = i.strip("\n")
            i = i.split("\t")

            gene = i[1]
            enh = i[3]
            if 25000 < float(i[4]) < 10000000:
                if gene not in enh_total.keys():
                    enh_total[gene] = [enh]
                else:
                    enh_total[gene].append(enh)

                if enh in align_score.keys():
                    if gene not in enh_conserv.keys():
                        enh_conserv[gene] = [float(align_score[enh])]
                    else:
                        enh_conserv[gene].append(float(align_score[enh]))

    output_file = path_evol + enh_name + "_" + data + "_summary_conserv.txt"
    output = open(output_file, 'w')
    if os.stat(output_file).st_size == 0:
        output.write("gene\tnb_total\tnb_seq_conserv\tnb_synt10M_conserv\tnb_synt2M_conserv\t"
                     "nb_contact_conserv\tmed_align_score\n")

    for gene in enh_total.keys():
        if gene in ortho.keys():
            nb_conserv = str(len(enh_conserv[gene])) if gene in enh_conserv.keys() else str(0)
            med_align = str(np.median(enh_conserv[gene])) if gene in enh_conserv.keys() else str(0)
            nb_synt10M = str(len(list(set(enh_synt10M[gene])))) if gene in enh_synt10M.keys() else str(0)
            nb_synt2M = str(len(list(set(enh_synt2M[gene])))) if gene in enh_synt2M.keys() else str(0)
            nb_contact = str(len(list(set(enh_cont[gene])))) if gene in enh_cont.keys() else str(0)

            output.write(gene + '\t' + str(len(set(enh_total[gene]))) + '\t' + nb_conserv + '\t' +
                         nb_synt10M + '\t' + nb_synt2M + '\t' + nb_contact + '\t' + med_align + '\n')

    output.close()


datas = ["original", "simulated"]
enh_data = ["CAGE", "ENCODE"]
if ref_sp == "human":
    enh_data.extend(["GRO_seq", "RoadMap"])

for dat in datas:
    for enh in enh_data:
        print("Running", ref_sp, "to", target_sp, "in", dat, "for", enh)
        summary(enh, dat)











