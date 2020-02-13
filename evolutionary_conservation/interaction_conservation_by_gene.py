#!/usr/bin/env python3
# coding=utf-8

import os
import numpy as np

# Conservation mouse interaction in human:
ref_sp = "human"
target_sp = "mouse"

path_data = "/home/laverre/Documents/Regulatory_Landscape/data/"
path_result = "/home/laverre/Documents/Regulatory_Landscape/result/"
seuil = 0.1


def running_all(data_ref, data_target):
    # Chromatin contact dictionary in both species
    def interaction_dict(sp, data):
        interaction = {}
        stats = {}
        if data == "_simul":
            input = "/Simulations/simulations_" + sp + "_10Mb_bin5kb_fragoverbin_chr_merged.txt_corrected2"
        else:
            input = "/all_interactions/all_interactions_chr_merged.txt_cell_names_corrected2"

        with open(path_data + sp + input) as f3:
            for i in f3.readlines()[1:]:
                i = i.strip("\n")
                i = i.split("\t")
                bait = (str(i[0]) + ":" + str(i[1]) + ":" + str(i[2]))
                PIR = (str(i[3]) + ":" + str(i[4]) + ":" + str(i[5]))

                if data == "":
                    nb_type = float(i[9])
                    median_strength = np.median(list(map(float, i[10].split(","))))
                else:
                    nb_type = median_strength = "NA"

                if i[0] == i[3]:
                    dist = float(i[7])
                    if 25000 < dist < 10000000:
                        if bait in interaction.keys():
                            interaction[bait].append(PIR)
                        else:
                            interaction[bait] = [PIR]

                else:
                    dist = "trans"

                stats[bait + '-' + PIR] = (str(dist), str(nb_type), str(median_strength), str(i[6]))

        return interaction, stats

    interaction_ref, stats_ref = interaction_dict(ref_sp, data_ref)
    interaction_target, stats_target = interaction_dict(target_sp, data_target)
    print("Interaction dictionary done !")

    # Genes to Bait
    def dict_gene2bait(sp, data):
        gene2bait = {}
        with open(path_result + "conservation/bait_composition_" + sp + data + "_merged.txt") as f3:
            for i in f3.readlines()[1:]:
                i = i.strip("\n")
                i = i.split("\t")
                bait = (str(i[0]) + ":" + str(i[1]) + ":" + str(i[2]))

                if i[len(i) - 1] != '':
                    genes = i[len(i) - 1].split(',')
                    for gene in genes:
                        if gene not in gene2bait.keys():
                            gene2bait[gene] = [bait]
                        elif bait not in gene2bait[gene]:
                            gene2bait[gene].append(bait)

        return gene2bait

    gene2bait_ref = dict_gene2bait(ref_sp, data_ref)
    gene2bait_target = dict_gene2bait(target_sp, data_target)
    print("Gene2Baits dictionary done !")

    # Orthologous gene one2one
    ortho = {}
    with open(path_data + "human2mouse_ortholog_one2one_genes_Ensembl94.txt") as f3:
        for i in f3.readlines()[1:]:
            i = i.strip("\n")
            i = i.split("\t")
            gene_origin = str(i[0]) if ref_sp == "human" else str(i[6])
            gene_target = str(i[6]) if ref_sp == "human" else str(i[0])

            ortho[gene_origin] = gene_target

    # Align score for each fragment of ref sp in target sp
    frag_conserv = {}
    if data_ref == "sample":
        file = path_result + "conservation/alignments/" + ref_sp + "/contacted_seq/AlignmentStatistics_Excluding_all_Exons_" + ref_sp + "2" + target_sp + "_samples_simulations.txt"
    else:
        file = path_result + "conservation/alignments/" + ref_sp + "/contacted_seq/AlignmentStatistics_Excluding_all_Exons_" + ref_sp + "2" + target_sp + "_merged_obs_simul_bait.txt"

    with open(file) as f1:
        for i in f1.readlines()[1:]:
            i = i.strip("\n")
            i = i.split("\t")
            origin_frag = i[0].strip(':+')
            origin_frag = origin_frag.strip(':-')
            frag_align = i[1].strip(':+')
            frag_align = frag_align.strip(':-')

            score_all_ungapped = int(i[4]) / (int(i[8]) + 1)

            if score_all_ungapped > seuil:
                frag_conserv[origin_frag] = (frag_align, score_all_ungapped)

    print("Score align : done ! ")

    # Overlap homologous fragment from ref sp to fragment in target sp
    overlap_target = {}
    with open(
            path_data + target_sp + "/overlap/" + ref_sp + "2" + target_sp + "_merged_overlap_" + target_sp + "_merged.txt") as f2:
        for i in f2.readlines()[1:]:
            i = i.strip("\n")
            i = i.split("\t")
            frag_align = str(i[1] + ':' + str(i[2]) + ':' + str(i[3]))
            if i[4] != "NA":
                overlap_target[frag_align] = i[4].split(',')  # Add all overlapping fragment

    print("Overlap aligned with target sp frag : done !")

    # Duplication score
    def dict_dupli(sp, data):
        dic = {}
        if data == "sample":
            file = path_result + "BLAT_duplication/" + sp + "_samples_simulations_BLAT_summary_0.8.txt"
        else:
            file = path_result + "BLAT_duplication/" + sp + "_merged_obs_simul_bait_BLAT_summary_0.8.txt"

        with open(file) as f3:
            for i in f3.readlines()[1:]:
                i = i.strip("\n")
                i = i.split("\t")
                frag = i[0]
                dic[frag] = int(i[4]) - 1

        return dic

    frag_dupli = dict_dupli(ref_sp, data_ref)
    frag_dupli_ortho = dict_dupli(target_sp, data_ref)  # data target !!!!!
    print("Score dupli : done ! ")

    print("Running conservation of interactions... ")
    conserv_contact = {}
    summary_conserv = {}

    print("Gene2bait ref:", len(gene2bait_ref), "Gene2bait target:", len(gene2bait_target))
    print("Bait ref:", len(interaction_ref.keys()), "Bait target:", len(interaction_target.keys()))

    for gene_ref in ortho.keys():
        nb_contact = nb_overlap = seq_conserv = seq_conserv_50 = contact_conserv = synt_conserv = 0
        # Are orthologous genes in baited fragments in both sides ?
        if gene_ref in gene2bait_ref.keys():
            if ortho[gene_ref] in gene2bait_target.keys():
                for bait_ref in gene2bait_ref[gene_ref]:
                    for bait_target in gene2bait_target[ortho[gene_ref]]:

                        # Are baited fragments in PC-HiC data in both sides ?
                        if bait_ref in interaction_ref.keys():
                            if bait_target in interaction_target.keys():

                                # Is contacted fragment conserved in target sp ?
                                for contact_ref in interaction_ref[bait_ref]:
                                    nb_contact += 1
                                    if contact_ref in frag_conserv.keys():
                                        seq_conserv += 1
                                        conserv_frag = frag_conserv[contact_ref][0]
                                        overlap_score = frag_conserv[contact_ref][1]
                                        if overlap_score > 0.5:
                                            seq_conserv_50 += 1

                                        # Is this fragment overlap with contacted fragment in target sp ?
                                        if conserv_frag in overlap_target.keys():
                                            for homolog_frag in overlap_target[conserv_frag]:
                                                nb_overlap += 1

                                                # Is this homologue fragment in synteny ?
                                                if bait_target.split(':')[0] == homolog_frag.split(':')[0]:  # same chromosome
                                                    midbait = (int(bait_target.split(':')[2]) + int(
                                                        bait_target.split(':')[1])) / 2
                                                    midfrag = (int(homolog_frag.split(':')[2]) + int(
                                                        homolog_frag.split(':')[1])) / 2
                                                    dist = abs(midbait - midfrag)
                                                    if dist < 10000000:
                                                        synt_conserv += 1

                                                # Is this homologue fragment in contact with orthologue gene ?
                                                if homolog_frag in interaction_target[bait_target]:
                                                    conserv_contact[gene_ref + '\t' + contact_ref] \
                                                        = (bait_ref, bait_target, homolog_frag, overlap_score,
                                                           "interaction_conserved")
                                                    contact_conserv += 1

                                                # It's in contact with an other bait/genes
                                                elif gene_ref + '\t' + contact_ref not in conserv_contact.keys():
                                                    conserv_contact[gene_ref + '\t' + contact_ref] \
                                                        = (
                                                    bait_ref, bait_target, homolog_frag, overlap_score, "not_good_gene")

                                                elif conserv_contact[gene_ref + '\t' + contact_ref][
                                                    4] != "interaction_conserved":
                                                    conserv_contact[gene_ref + '\t' + contact_ref] \
                                                        = (
                                                    bait_ref, bait_target, homolog_frag, overlap_score, "not_good_gene")

                                        # Conserved but not in contact
                                        elif gene_ref + '\t' + contact_ref not in conserv_contact.keys():
                                            conserv_contact[gene_ref + '\t' + contact_ref] = (
                                            bait_ref, bait_target, conserv_frag, overlap_score, "not_in_contact")
                                        elif conserv_contact[gene_ref + '\t' + contact_ref][
                                            4] != "interaction_conserved" and \
                                                conserv_contact[gene_ref + '\t' + contact_ref][4] != "not_good_gene":
                                            conserv_contact[gene_ref + '\t' + contact_ref] = (
                                            bait_ref, bait_target, conserv_frag, overlap_score, "not_in_contact")

                                    # Not conserved
                                    elif gene_ref + '\t' + contact_ref not in conserv_contact.keys():
                                        conserv_contact[gene_ref + '\t' + contact_ref] = (
                                        bait_ref, bait_target, "NA", "NA", "not_conserved")

                                summary_conserv[gene_ref] = [nb_contact, nb_overlap, seq_conserv, seq_conserv_50,
                                                             synt_conserv, contact_conserv]

    output_file = "/conservation/interaction_conservation/" + ref_sp + data_ref + "_merged_to_" + target_sp + data_target + "_by_gene.txt2"

    output = open(path_result + output_file, 'w')
    if os.stat(path_result + output_file).st_size == 0:
        output.write(
            "gene\tbait\tcontact\tdist\tnb_cell\tmedian_strength\tbaited_contact\tduplication\tgene_ortho\tbait_ortho\tcontact_ortho"
            "\tdist_ortho\tnb_cell_ortho\tmedian_strength_ortho\tbaited_contact_ortho\tduplication_ortho\toverlap\tstatus\tconserved_frag\n")

    for inter in conserv_contact.keys():
        gene_ref = inter.split('\t')[0]
        contact_ref = inter.split('\t')[1]
        bait_target = conserv_contact[inter][1]
        conserv_frag = conserv_contact[inter][2]

        if contact_ref not in frag_dupli.keys():  # not find by BLAT
            frag_dupli[contact_ref] = 0

        if conserv_frag not in frag_dupli_ortho.keys():  # not find by BLAT
            frag_dupli_ortho[conserv_contact[inter][2]] = 0

        if contact_ref not in frag_conserv.keys():  # not find by liftOver
            frag_conserv[contact_ref] = ("NA", 0)

        interaction_ref = conserv_contact[inter][0] + '-' + contact_ref
        interaction_target = bait_target + '-' + conserv_frag

        if interaction_target not in stats_target.keys():
            if bait_target[0] == conserv_frag.split(':')[0]:  # same chromosome
                midbait = (int(bait_target[2]) + int(bait_target[1])) / 2
                midfrag = (int(conserv_frag.split(':')[2]) + int(conserv_frag.split(':')[1])) / 2
                dist = abs(midbait - midfrag)
                stats_target[interaction_target] = (str(dist), "NA", "NA", "NA")

            else:
                stats_target[interaction_target] = ("trans", "NA", "NA", "NA")

        output.write(gene_ref + '\t' + str(conserv_contact[inter][0]) + '\t' + contact_ref + '\t' +
                     str("\t".join(stats_ref[interaction_ref])) + '\t' + str(frag_dupli[contact_ref]) + '\t' +
                     str(ortho[gene_ref]) + '\t' + str(conserv_contact[inter][1]) + '\t' +
                     str(conserv_contact[inter][2]) + '\t' + str("\t".join(stats_target[interaction_target])) + '\t' +
                     str(frag_dupli_ortho[conserv_contact[inter][2]]) + '\t' + str(conserv_contact[inter][3]) + '\t' +
                     str(conserv_contact[inter][4]) + '\t' + str(frag_conserv[contact_ref][0]) + '\n')

    output.close()

    output_file = "/conservation/" + ref_sp + data_ref + "_merged_to_" + target_sp + data_target + "_summary_by_gene.txt"
    output_summary = open(path_result + output_file, 'w')
    if os.stat(path_result + output_file).st_size == 0:
        output_summary.write(
            "gene\tnb_contact\toverlap_target\tseq_conserv\tseq_conserv_50\tsynt_conserv\tcontact_conserv\n")

    for gene in summary_conserv.keys():
        output_summary.write(gene + '\t' + str("\t".join(str(x) for x in summary_conserv[gene])) + '\n')

    output_summary.close()


data1 = ["", "_simul"]  # or "" for observed data
data2 = ["", "_simul"]

for ref in data1:
    for target in data2:
        print("Origin sp:", ref_sp, "; data:", ref, "; to data :", target)
        running_all(ref, target)
