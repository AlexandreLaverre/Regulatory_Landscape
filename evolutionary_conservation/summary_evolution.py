#!/usr/bin/env python3
# coding=utf-8

import os
import numpy as np

ref_sp = "human"
target_sp = "mouse" if ref_sp == "human" else "human"
treshold = 0.4
min_dist = 25000

pathData = "/home/laverre/Manuscript/SupplementaryDataset1/" + ref_sp + "/"
pathSimulation = "/home/laverre/Manuscript/SupplementaryDataset2/" + ref_sp + "/"
pathEnhContact = "/home/laverre/Manuscript/SupplementaryDataset4/" + ref_sp + "/"

pathEvolution = "/home/laverre/Manuscript/SupplementaryDataset7/" + ref_sp
pathAlignment = pathEvolution + "/sequence_conservation/enhancers/"
pathSynteny = pathEvolution + "/synteny_conservation/"
pathContact = pathEvolution + "/contact_conservation/"
pathOutput = pathEvolution + "/evolution_summary_by_gene/"

path_annot = "/home/laverre/Data/Regulatory_landscape/result/Supplementary_dataset3_annotations/" + ref_sp + "/"

########################################### Orthologous gene one2one ###############################################
ortho = {}
with open(pathEvolution + "/gene_orthology/" + ref_sp + "2" + target_sp + "_orthologue.txt") as f3:
    for i in f3.readlines()[1:]:
        i = i.strip("\n")
        i = i.split("\t")
        if i[9] == "ortholog_one2one":
            ortho[str(i[0])] = str(i[5])

#################################################### Bait to Genes ####################################################
bait2gene = {}
with open(path_annot + "/restriction_fragments/bait_overlap_TSS_1kb.txt") as f1:
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
        infile = pathData + "/all_interactions.txt"
    else:
        infile = pathSimulation + "/simulated_all_interactions.txt"

    coverage = {}
    gene_contact = {}
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
                    else [colnames.index("ESC"), colnames.index("ESC_18"), colnames.index("ESC_wild")]
            elif sample == "Bcell":
                ref = [colnames.index("TB"), colnames.index("NB")] if ref_sp == "human" \
                    else [colnames.index("preB_aged"), colnames.index("preB_young")]
            else:
                ref = [1]

            # Get only specific contacts
            not_ref = [x for x in range(8, len(colnames)) if x not in ref]
            #and all(i[y] == "NA" for y in not_ref)

            if any(i[x] != "NA" for x in ref):
                if i[0] == i[3]:    # cis-interaction
                    if min_dist <= float(i[7]) <= 10000000:
                        if i[6] == "unbaited":

                            if bait in bait2gene.keys():
                                for gene in bait2gene[bait]:
                                    if gene not in gene_contact.keys():
                                        gene_contact[gene] = [contacted]
                                        coverage[gene] = [contacted_size]
                                    elif contacted not in gene_contact[gene]:
                                        gene_contact[gene].append(contacted)
                                        coverage[gene].append(contacted_size)

    for gene in coverage.keys():
        coverage[gene] = str(sum(size for size in coverage[gene]))

    return coverage


############################################# Enhancers alignments ###################################################
def enh_score(enh_name):
    duplication = {}
    align = {}
    with open(path_annot + "/" + enh_name + "/" + enh_name + "_BLAT_summary_0.8.txt") as f1:
        for i in f1.readlines()[1:]:
            i = i.strip("\n")
            i = i.split("\t")
            enh = i[0]

            if int(i[4]) == 1:
                duplication[enh] = i[4]

    with open(pathAlignment + enh_name + "/AlignmentStatistics_Excluding_Exons_" + ref_sp + "2" + target_sp + ".txt") as f1:
        for i in f1.readlines()[1:]:
            i = i.strip("\n")
            i = i.split("\t")
            enh = i[0].strip(":+")

            try:
                align_score = int(i[6]) / int(i[9])  # FilteredUngappedLength / FilteredAlignmentLength
            except ZeroDivisionError:
                align_score = 0

            if enh in duplication.keys(): #align_score > treshold and
                align[enh] = align_score

    overlap_target = {}
    if enh_name in ["FANTOM5", "ENCODE"]:
        with open(pathAlignment + enh_name + "/enhancer_overlap_target_enhancer.txt") as f1:
            for i in f1.readlines()[1:]:
                i = i.strip("\n")
                i = i.split("\t")

                enh = i[0]
                if enh in duplication.keys() and i[4] != "NA":
                    overlap_target[enh] = i[4]

    return duplication, align, overlap_target


############################################# Synteny conservation ###################################################
def enh_synteny(enh_name, data, enh_conserved):
    synteny_10M = {}
    synteny_2M = {}
    with open(pathSynteny + enh_name + "/" + ref_sp + "2" + target_sp + "_" + data + "_synteny.txt") as f1:
        for i in f1.readlines()[1:]:
            i = i.strip("\n")
            i = i.split("\t")
            gene = i[0]
            enh = i[2]
            gene_enh_dist = float(i[3])

            if gene_enh_dist > min_dist:
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
    with open(pathContact + enh_name + "/" + ref_sp + "_" + data + "2" + target_sp + "_" + data + ".txt") as f1:

        colnames = f1.readline().strip("\n")
        colnames = colnames.split("\t")

        for i in f1.readlines():
            i = i.strip("\n")
            i = i.split("\t")
            gene = i[0]
            enh = i[1]
            dist = float(i[2])

            if dist > min_dist:
                if sample == "pre_adipo":
                    ref = [colnames.index("pre_adipo")]
                    target = [colnames.index("preadip_D0"), colnames.index("preadip_D2"), colnames.index("preadip_4H")]
                elif sample == "ESC":
                    ref = [colnames.index("hESC")]
                    target = [colnames.index("ESC"), colnames.index("ESC_18"),  colnames.index("ESC_wild")]
                elif sample == "Bcell":
                    ref = [colnames.index("TB"), colnames.index("NB")]
                    target = [colnames.index("preB_aged"), colnames.index("preB_young")]
                else:
                    target = [1]
                    if int(i[4]) > 1:
                        ref = [1]
                    else:
                        continue

                col_ref = range(colnames.index("Bcell"), colnames.index("NCD8")+1)
                col_target = range(colnames.index("EpiSC"), colnames.index("preadip_4H")+1)

                if any(i[x] != '0' for x in ref) and any(i[x] != '0' for x in target):
                    if enh in duplication.keys():
                        if enh in align_score.keys() and align_score[enh] > treshold:
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
    coverage = cover(data, sample)
    duplication, align_score, overlap_target = enh_score(enh_name)
    enh_cont, contact_overlap = enh_contact(enh_name, data, sample, duplication, align_score, overlap_target)
    enh_total = {}
    dist = {}
    enh_total_length = {}
    enh_conserv = {}
    enh_conserv_list = {}
    with open(pathEnhContact + "/" + enh_name + "/gene_enhancer_contacts_" + data + "_interactions.txt") as f1:
        colnames = f1.readline().strip("\n")
        colnames = colnames.split("\t")

        for i in f1.readlines():
            i = i.strip("\n")
            i = i.split("\t")
            gene = i[0]
            enh = i[1]
            gene_enh_dist = float(i[4])
            enh_length = int(enh.split(":")[2]) - int(enh.split(":")[1])

            if gene_enh_dist > min_dist:
                if gene not in dist.keys():
                    dist[gene] = [float(i[4])]
                else:
                    dist[gene].append(float(i[4]))

                if sample == "pre_adipo":
                    ref = [colnames.index("pre_adipo")] if ref_sp == "human" \
                        else [colnames.index("preadip_D0"), colnames.index("preadip_D2"), colnames.index("preadip_4H")]
                elif sample == "ESC":
                    ref = [colnames.index("hESC")] if ref_sp == "human" \
                        else [colnames.index("ESC"), colnames.index("ESC_18"), colnames.index("ESC_wild")]
                elif sample == "Bcell":
                    ref = [colnames.index("TB"), colnames.index("NB")] if ref_sp == "human" \
                        else [colnames.index("preB_aged"), colnames.index("preB_young")]
                else:
                    if int(i[6]) > 1:
                        ref = [1]
                    else:
                        continue

                if any(i[x] != "nan" for x in ref):
                    if enh in duplication.keys():
                        if enh not in align_score.keys():
                            align_score[enh] = 0

                        if gene not in enh_total.keys():
                            enh_total[gene] = [float(align_score[enh])]
                            enh_total_length[gene] = [enh_length]
                        elif enh not in enh_total[gene]:
                            enh_total[gene].append(float(align_score[enh]))
                            enh_total_length[gene].append(enh_length)

                        if enh in align_score.keys() and align_score[enh] > treshold:
                            if gene not in enh_conserv.keys():
                                enh_conserv[gene] = [float(align_score[enh])]
                                enh_conserv_list[gene] = [enh]
                            elif enh not in enh_conserv_list[gene]:
                                enh_conserv[gene].append(float(align_score[enh]))
                                enh_conserv_list[gene].append(enh)

    enh_synt10M, enh_synt2M = enh_synteny(enh_name, data, enh_conserv_list)

    output_file = pathOutput + enh_name + "/" + data + "_evolution_summary_" + sample + "_" + str(treshold) + ".txt"
    output = open(output_file, 'w')
    if os.stat(output_file).st_size == 0:
        output.write("gene\tnb_total\tnb_seq_conserv\tnb_synt10M_conserv\tnb_synt2M_conserv\t"
                     "nb_contact_conserv\tmed_align_score\tmean_align_score\tcontacted_coverage\ttotal_enh_length\tmedian_dist")

        if enh_name in ["FANTOM5", "ENCODE"]:
            output.write("\tnb_overlap_target\n")
        else:
            output.write("\n")

    for gene in coverage.keys():
        if gene in ortho.keys() and gene in enh_total.keys() and len(enh_total[gene]) > 1:

            nb_total = str(len(enh_total[gene])) if gene in enh_total.keys() else str(0)
            enh_cover = str(sum(length for length in enh_total_length[gene])) if gene in enh_total_length.keys() else str(0)
            nb_conserv = str(len(enh_conserv[gene])) if gene in enh_conserv.keys() else str(0)
            med_align = str(np.median(enh_total[gene])) if gene in enh_total.keys() else str(0)
            mean_align = str(np.mean(enh_total[gene])) if gene in enh_total.keys() else str(0)
            nb_synt10M = str(len(enh_synt10M[gene])) if gene in enh_synt10M.keys() else str(0)
            nb_synt2M = str(len(enh_synt2M[gene])) if gene in enh_synt2M.keys() else str(0)
            nb_contact = str(len(enh_cont[gene])) if gene in enh_cont.keys() else str(0)
            median_dist = str(np.median(dist[gene])) if gene in dist.keys() else str(0)

            output.write(gene + '\t' + nb_total + '\t' + nb_conserv + '\t' + nb_synt10M + '\t' + nb_synt2M + '\t' +
                         nb_contact + '\t' + med_align + '\t' + mean_align + '\t' + coverage[gene] + '\t' + enh_cover + '\t' + median_dist)

            if enh_name in ["FANTOM5", "ENCODE"]:
                nb_overlap_target = str(len(contact_overlap[gene])) if gene in contact_overlap.keys() else str(0)
                output.write("\t" + nb_overlap_target + "\n")

            else:
                output.write("\n")

    output.close()


datas = ["original", "simulated"]
enh_data = ["FANTOM5", "ENCODE"]
samples = ["all", "Bcell", "pre_adipo", "ESC"]
if ref_sp == "human":
     enh_data.extend(["FOCS_GRO_seq", "RoadmapEpigenomics"])

for enh_dat in enh_data:
    for dat in datas:
        for samp in samples:
            print("Running", ref_sp, "to", target_sp, "in", dat, "for", enh_dat, "in", samp, "sample")
            summary(enh_dat, dat, samp)

