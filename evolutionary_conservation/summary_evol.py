#!/usr/bin/env python3
# coding=utf-8

import os
import numpy as np

ref_sp = "human"
target_sp = "mouse" if ref_sp == "human" else "human"
seuil = 0.6

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

#################################################### Bait to Genes ####################################################
bait2gene = {}
with open(path_annot + "bait_overlap_TSS_1kb.txt") as f1:
    for i in f1.readlines()[1:]:
        i = i.strip("\n")
        i = i.split("\t")
        bait = (str(i[1]) + ":" + str(i[2]) + ":" + str(i[3]))

        if i[5] != "NA":
            genes = i[5].split(',')
            for gene in genes:
                if bait not in bait2gene.keys():
                    bait2gene[bait] = [gene]
                else:
                    bait2gene[bait].append(gene)


######################################## Coverage of contacted region per gene ########################################
def cover(data, sample):
    if data == "original":
        infile = path + "/Supplementary_dataset1_original_interactions/" + ref_sp + "/all_interactions.txt"
    else:
        infile = path + "/Supplementary_dataset2_simulated_interactions/" + ref_sp + "/simulated_all_interactions.txt"

    coverage = {}
    gene_contact = {}
    dist = {}
    with open(infile) as f1:
        colnames = f1.readline().strip("\n")
        colnames = colnames.split("\t")

        for i in f1.readlines():
            i = i.strip("\n")
            i = i.split("\t")

            bait = (str(i[0]) + ":" + str(i[1]) + ":" + str(i[2]))
            contacted = (str(i[3]) + ":" + str(i[4]) + ":" + str(i[5]))
            contacted_size = int(i[5]) - int(i[4])

            if sample == "pre_adipo":
                ref = [colnames.index("pre_adipo")] if ref_sp == "human" \
                    else [colnames.index("preadip_D0"), colnames.index("preadip_D2"), colnames.index("preadip_4H")]
            elif sample == "ESC":
                ref = [colnames.index("hESC")] if ref_sp == "human" \
                    else [colnames.index("ESC"), colnames.index("ESC_18"), colnames.index("ESC_NKO"),  colnames.index("ESC_wild")]
            elif sample == "Bcell":
                ref = [colnames.index("Bcell"), colnames.index("TB")] if ref_sp == "human" \
                    else [colnames.index("preB_aged"), colnames.index("preB_young")]
            else:
                ref = [1]

            if any(i[x] != "NA" for x in ref):
                if i[0] == i[3]:    # cis-interaction
                    if 25000 <= float(i[7]) <= 10000000:
                        if i[6] == "unbaited":

                            if bait in bait2gene.keys():
                                for gene in bait2gene[bait]:
                                    if gene not in gene_contact.keys():
                                        gene_contact[gene] = [contacted]
                                        coverage[gene] = [contacted_size]
                                        dist[gene] = [float(i[7])]
                                    elif contacted not in gene_contact[gene]:
                                        gene_contact[gene].append(contacted)
                                        coverage[gene].append(contacted_size)
                                        dist[gene].append(float(i[7]))

    for gene in coverage.keys():
        coverage[gene] = str(sum(size for size in coverage[gene]))

    return coverage, dist


############################################# Enhancers alignments ###################################################
def enh_score(enh_name):
    duplication = {}
    align = {}
    with open(path_annot + enh_name + "_BLAT_summary_0.8.txt") as f1:
        for i in f1.readlines()[1:]:
            i = i.strip("\n")
            i = i.split("\t")
            enh = i[0]

            if int(i[4]) == 1:
                duplication[enh] = i[4]

    with open(path_evol + "sequence_conservation/" + enh_name +
              "/AlignmentStatistics_Excluding_all_Exons_" + ref_sp + "2" + target_sp + "_" + enh_name + ".txt") as f1:
        for i in f1.readlines()[1:]:
            i = i.strip("\n")
            i = i.split("\t")
            enh = i[0].strip(":+")

            try:
                align_score = int(i[6]) / int(i[9])  # FilteredUngappedLength / FilteredAlignmentLength
            except ZeroDivisionError:
                align_score = 0

            if align_score > seuil and enh in duplication.keys():
                align[enh] = align_score

    overlap_target = {}
    if enh_name in ["CAGE", "ENCODE"]:
        with open(path_evol + "sequence_conservation/" + enh_name + "/" +
                  enh_name + "_lifted_overlap_" + enh_name + "_target.txt") as f1:
            for i in f1.readlines()[1:]:
                i = i.strip("\n")
                i = i.split("\t")

                enh = i[0]
                if enh in duplication.keys() and i[4] != "NA":
                    overlap_target[enh] = i[4]

    return duplication, align, overlap_target

    # with open(path_evol + "/sequence_conservation/PhastCons/PhastCons_vertebrates_" +
    #           enh_name + "_MaskedExons_Ensembl94.txt") as f1:
    #     for i in f1.readlines()[1:]:
    #         i = i.strip("\n")
    #         i = i.split("\t")
    #
    #         enh_origin = i[0].strip(":+")
    #         if i[4] == "NA":
    #             align_score = 0
    #         else:
    #             align_score = float(i[4])
    #
    #         if align_score > seuil:
    #             align[enh_origin] = align_score
    #
    # return align


############################################# Synteny conservation ###################################################
def enh_synteny(enh_name, data, enh_conserved):
    synteny_10M = {}
    synteny_2M = {}
    with open(path_evol + "synteny_conservation/" + enh_name + "/" + ref_sp + "2" + target_sp + "_"
              + enh_name + "_" + data + "_synteny.txt") as f1:
        for i in f1.readlines()[1:]:
            i = i.strip("\n")
            i = i.split("\t")
            gene = i[0]
            enh = i[2]

            if gene in enh_conserved.keys() and enh in enh_conserved[gene]:
                if i[7] != "trans" and float(i[7]) < 10000000:          # cis interaction
                    if gene not in synteny_10M.keys():
                        synteny_10M[gene] = [enh]
                    elif enh not in synteny_10M[gene]:
                        synteny_10M[gene].append(enh)

                    if float(i[7]) < 2000000:
                        if gene not in synteny_2M.keys():
                            synteny_2M[gene] = [enh]
                        elif enh not in synteny_2M[gene]:
                            synteny_2M[gene].append(enh)

    return synteny_10M, synteny_2M


############################################# Contact conservation ###################################################
def enh_contact(enh_name, data, sample, duplication, align_score, overlap_target):
    contact = {}
    contact_overlap = {}
    with open(path_evol + "contact_conservation/" + enh_name + "/" + ref_sp + "2" + target_sp +
              "_" + data + ".txt") as f1:

        colnames = f1.readline().strip("\n")
        colnames = colnames.split("\t")

        for i in f1.readlines():
            i = i.strip("\n")
            i = i.split("\t")
            gene = i[0]
            enh = i[1]

            if sample == "pre_adipo":
                ref = [colnames.index("pre_adipo")]
                target = [colnames.index("preadip_D0"), colnames.index("preadip_D2"), colnames.index("preadip_4H")]
            elif sample == "ESC":
                ref = [colnames.index("hESC")]
                target = [colnames.index("ESC"), colnames.index("ESC_18"), colnames.index("ESC_NKO"),  colnames.index("ESC_wild")]
            elif sample == "Bcell":
                ref = [colnames.index("Bcell"), colnames.index("TB")]
                target = [colnames.index("preB_aged"), colnames.index("preB_young")]
            else:
                ref = [1]
                target = [1]

            if any(i[x] != '0' for x in ref) and any(i[x] != '0' for x in target):
                if enh in duplication.keys():
                    if enh in align_score.keys():
                        # Contact dict
                        if gene not in contact.keys():
                            contact[gene] = [enh]
                        elif enh not in contact[gene]:
                            contact[gene].append(enh)

                        if enh_name in ["CAGE", "ENCODE"]:
                            # Contact and overlap dict
                            if enh in overlap_target.keys():
                                if gene not in contact_overlap.keys():
                                    contact_overlap[gene] = [enh]
                                elif enh not in contact_overlap[gene]:
                                    contact_overlap[gene].append(enh)

    return contact, contact_overlap


############################################# Nb contact ###################################################
def summary(enh_name, data, sample):
    coverage, dist = cover(data, sample)
    duplication, align_score, overlap_target = enh_score(enh_name)
    enh_cont, contact_overlap = enh_contact(enh_name, data, sample, duplication, align_score, overlap_target)
    enh_total = {}
    enh_total_length = {}
    enh_conserv = {}
    enh_conserv_list = {}
    with open(path_contact + "gene_" + enh_name + "_enhancers_" + data + "_interactions.txt") as f1:
        colnames = f1.readline().strip("\n")
        colnames = colnames.split("\t")

        for i in f1.readlines():
            i = i.strip("\n")
            i = i.split("\t")
            gene = i[1]
            enh = i[3]
            enh_length = int(enh.split(":")[2]) - int(enh.split(":")[1])

            if sample == "pre_adipo":
                ref = [colnames.index("pre_adipo")] if ref_sp == "human" \
                    else [colnames.index("preadip_D0"), colnames.index("preadip_D2"), colnames.index("preadip_4H")]
            elif sample == "ESC":
                ref = [colnames.index("hESC")] if ref_sp == "human" \
                    else [colnames.index("ESC"), colnames.index("ESC_18"), colnames.index("ESC_NKO"),  colnames.index("ESC_wild")]
            elif sample == "Bcell":
                ref = [colnames.index("Bcell"), colnames.index("TB")] if ref_sp == "human" \
                    else [colnames.index("preB_aged"), colnames.index("preB_young")]
            else:
                ref = [1]

            if any(i[x] != "nan" for x in ref):
                if enh in duplication.keys():
                    if gene not in enh_total.keys():
                        enh_total[gene] = [enh]
                        enh_total_length[gene] = [enh_length]
                    elif enh not in enh_total[gene]:
                        enh_total[gene].append(enh)
                        enh_total_length[gene].append(enh_length)

                    if enh in align_score.keys():
                        if gene not in enh_conserv.keys():
                            enh_conserv[gene] = [float(align_score[enh])]
                            enh_conserv_list[gene] = [enh]
                        elif enh not in enh_conserv_list[gene]:
                            enh_conserv[gene].append(float(align_score[enh]))
                            enh_conserv_list[gene].append(enh)

    enh_synt10M, enh_synt2M = enh_synteny(enh_name, data, enh_conserv_list)

    output_file = path_evol + enh_name + "_" + data + "_summary_conserv_" + sample + "_" + str(seuil) + ".txt"
    output = open(output_file, 'w')
    if os.stat(output_file).st_size == 0:
        output.write("gene\tnb_total\tnb_seq_conserv\tnb_synt10M_conserv\tnb_synt2M_conserv\t"
                     "nb_contact_conserv\tmed_align_score\tcontacted_coverage\ttotal_enh_length\tmedian_dist")

        if enh_name in ["CAGE", "ENCODE"]:
            output.write("\tnb_overlap_target\n")
        else:
            output.write("\n")

    for gene in coverage.keys():
        if gene in ortho.keys():
            nb_total = str(len(enh_total[gene])) if gene in enh_total.keys() else str(0)
            enh_cover = str(sum(length for length in enh_total_length[gene])) if gene in enh_total_length.keys() else str(0)
            nb_conserv = str(len(enh_conserv[gene])) if gene in enh_conserv.keys() else str(0)
            med_align = str(np.median(enh_conserv[gene])) if gene in enh_conserv.keys() else str(0)
            nb_synt10M = str(len(enh_synt10M[gene])) if gene in enh_synt10M.keys() else str(0)
            nb_synt2M = str(len(enh_synt2M[gene])) if gene in enh_synt2M.keys() else str(0)
            nb_contact = str(len(enh_cont[gene])) if gene in enh_cont.keys() else str(0)
            median_dist = str(np.median(dist[gene]))

            output.write(gene + '\t' + nb_total + '\t' + nb_conserv + '\t' + nb_synt10M + '\t' + nb_synt2M + '\t' +
                         nb_contact + '\t' + med_align + '\t' + coverage[gene] + '\t' + enh_cover + '\t' + median_dist)

            if enh_name in ["CAGE", "ENCODE"]:
                nb_overlap_target = str(len(contact_overlap[gene])) if gene in contact_overlap.keys() else str(0)
                output.write("\t" + nb_overlap_target + "\n")

            else:
                output.write("\n")

    output.close()


datas = ["original", "simulated"]
enh_data = ["CAGE", "ENCODE"]
samples = ["all"] #, "pre_adipo", "Bcell", "ESC"]
if ref_sp == "human":
     enh_data.extend(["GRO_seq", "RoadMap"])

for enh_dat in enh_data:
    for dat in datas:
        for samp in samples:
            print("Running", ref_sp, "to", target_sp, "in", dat, "for", enh_dat, "in", samp, "sample")
            summary(enh_dat, dat, samp)

