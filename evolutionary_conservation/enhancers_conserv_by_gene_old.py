#!/usr/bin/env python3
# coding=utf-8

import os
import numpy as np
import itertools

# Conservation mouse interaction in human:
ref_sp = "mouse"
target_sp = "human"

path_data = "/home/laverre/Documents/Regulatory_Landscape/data/"
path_result = "/home/laverre/Documents/Regulatory_Landscape/result/"
path_alignment = path_result + "conservation/alignments/" + ref_sp + "/"
seuil = 0.1


def running_all(data_ref, data_target):

    ################################### Chromatin contact dictionary in both species ###################################
    def interaction_dict(sp, data):
        interaction = {}
        stats = {}
        if data == "_simul":
            input = "/Simulations/simulations_" + sp + "_10Mb_bin5kb_fragoverbin_chr_merged.txt_corrected2"
        elif data == "sample":
            input = "/Simulations/samples_simulated_all_interactions_bait_all_merged.txt_corrected2"
        else:
            input = "/all_interactions/all_interactions_chr_merged.txt_cell_names_corrected2"

        with open(path_data + sp + input) as f3:
            for i in f3.readlines()[1:]:
                i = i.strip("\n")
                i = i.split("\t")
                bait = (str(i[0]) + ":" + str(i[1]) + ":" + str(i[2]))
                PIR = (str(i[3]) + ":" + str(i[4]) + ":" + str(i[5]))
                length = (int(i[5]) - int(i[4]))

                if data == "":
                    nb_type = float(i[9])
                    median_strength = np.median(list(map(float, i[10].split(","))))
                else:
                    nb_type = median_strength = "NA"

                if i[0] == i[3]:
                    dist = float(i[7])
                    if 25000 < dist < 10000000:
                        if sp == ref_sp:
                            if bait in interaction.keys():
                                interaction[bait].append((PIR, length))
                            else:
                                interaction[bait] = [(PIR, length)]
                        else:
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

    ############################################### Genes to Bait ######################################################
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

    ########################################### Orthologous gene one2one ###############################################
    ortho = {}
    with open(path_data + "human2mouse_ortholog_one2one_genes_Ensembl94.txt") as f3:
        for i in f3.readlines()[1:]:
            i = i.strip("\n")
            i = i.split("\t")
            gene_origin = str(i[0]) if ref_sp == "human" else str(i[6])
            gene_target = str(i[6]) if ref_sp == "human" else str(i[0])

            ortho[gene_origin] = gene_target

    ######################################## Align score for each enhancer ############################################
    def enhancers_conserv(enhancer):
        enh_conserv = {}
        file = path_alignment + enhancer + "/AlignmentStatistics_Excluding_all_Exons_" + ref_sp + "2" + target_sp + "_" + enhancer + ".txt"
        with open(file) as f1:
            for i in f1.readlines()[1:]:
                i = i.strip("\n")
                i = i.split("\t")
                origin_enh = i[0].strip(':+')
                origin_enh = origin_enh.strip(':-')
                enh_align = i[1].strip(':+')
                enh_align = enh_align.strip(':-')

                score_filtr_ungapped = int(i[6]) / (int(i[9]) + 1)  # Total = int(i[4]) / (int(i[8]) + 1)
                if score_filtr_ungapped > seuil:
                    enh_conserv[origin_enh] = (enh_align, score_filtr_ungapped)

            return(enh_conserv)

    CAGE_conserv = enhancers_conserv("CAGE")
    
    if ref_sp == "human":
        ENCODE_conserv = enhancers_conserv("ENCODE")
        RoadMap_conserv = enhancers_conserv("RoadMap")
        GRO_seq_conserv = enhancers_conserv("GRO_seq")

    print("Score align : done ! ")

    ######################################### Overlap enhancers and fragments ##########################################
    def enhancers(sp, file):
        enh = {}
        if sp == ref_sp:
            overlap = path_data + ref_sp + "/overlap/" + ref_sp + file
        else:
            overlap = path_data + target_sp + "/overlap/" + target_sp + file
        with open(overlap) as f2:
            for i in f2.readlines()[1:]:
                i = i.strip("\n")
                i = i.split("\t")
                frag = (str(i[1]) + ':' + str(i[2]) + ':' + str(i[3]))
                if sp == ref_sp:
                    if i[4] != "NA":
                        enh[frag] = i[4].split(',')  # Add all enhancers
                else:
                    if i[4] != "NA":
                        for enh_lift in i[4].split(','):
                            if enh_lift not in enh.keys():
                                enh[enh_lift] = [frag]
                            else:
                                enh[enh_lift].append(frag)  # Add all overlapped fragments
        return enh

    CAGE_ref = enhancers(ref_sp, "_merged_overlap_CAGE.txt")
    CAGE_target = enhancers(target_sp, str("_merged_overlap_"+ref_sp+"2"+target_sp+"_CAGE.txt"))

    dic_enh = {"CAGE": [CAGE_conserv, CAGE_ref, CAGE_target]} #, "ENCODE": [ENCODE_conserv, ENCODE_ref, ENCODE_target]}

    if ref_sp == "human":
        ENCODE_ref = enhancers(ref_sp, "_merged_overlap_ENCODE.txt")
        ENCODE_target = enhancers(target_sp, str("_merged_overlap_" + ref_sp + "2" + target_sp + "_ENCODE.txt"))
        RoadMap_ref = enhancers(ref_sp, "_merged_overlap_RoadMap.txt")
        RoadMap_target = enhancers(target_sp, str("_merged_overlap_"+ref_sp+"2"+target_sp+"_RoadMap.txt"))
        GRO_seq_ref = enhancers(ref_sp, "_merged_overlap_GRO_seq.txt")
        GRO_seq_target = enhancers(target_sp, str("_merged_overlap_"+ref_sp+"2"+target_sp+"_GRO_seq.txt"))
        dic_enh["RoadMap"] = [RoadMap_conserv, RoadMap_ref, RoadMap_target]
        dic_enh["GRO_seq"] = [GRO_seq_conserv, GRO_seq_ref, GRO_seq_target]

    print("Overlap aligned with target sp frag : done !")
    ################################################ Duplication score #################################################
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

    ########################################## Conservation by enhancer ################################################
    def conserv_enh(enh_data, origin_contact, origin_gene):
        nb_enh = nb_enh_conserv = nb_enh_overlap = nb_enh_synt = nb_enh_contact = 0

        dic_enh_conserv = dic_enh[enh_data][0]
        dic_enh_ref = dic_enh[enh_data][1]
        dic_enh_target = dic_enh[enh_data][2]

        if origin_contact in dic_enh_ref.keys():
            for enh in dic_enh_ref[origin_contact]:
                synt = contact = "F"
                nb_enh += 1

                # Does this enhancer conserved ?
                if enh in dic_enh_conserv.keys():
                    nb_enh_conserv += 1
                    enh_lift = dic_enh_conserv[enh][0]

                    # Does this homologue enhancer conserved in synteny ?
                    for bait_target in gene2bait_target[ortho[origin_gene]]:
                        if bait_target.split(':')[0] == enh_lift.split(':')[0]:  # same chromosome
                            midbait = (int(bait_target.split(':')[2]) + int(bait_target.split(':')[1])) / 2
                            midfrag = (int(enh_lift.split(':')[2]) + int(enh_lift.split(':')[1])) / 2
                            dist = abs(midbait - midfrag)
                            if dist < 10000000:
                                synt = "T"

                    # Is this lifted enhancer overlap with contacted fragment in target sp ?
                    if enh_lift in dic_enh_target.keys():
                        for contact_target in dic_enh_target[enh_lift]:
                            nb_enh_overlap += 1

                            # Does target bait in PC-HiC data ?
                            for bait_target in gene2bait_target[ortho[gene_ref]]:
                                if bait_target in interaction_target.keys():
                                    # Is this homologue enh in contact with orthologue gene ?
                                    if contact_target in interaction_target[bait_target]:
                                        contact = "T"

                if synt == "T": nb_enh_synt += 1
                if contact == "T": nb_enh_contact += 1

        return [nb_enh, nb_enh_conserv, nb_enh_overlap, nb_enh_synt, nb_enh_contact]

    ############################################# Conservation analyses ################################################
    print("Running conservation of interactions... ")
    conserv_contact = {}
    summary_conserv = {}

    print("Gene2bait ref:", len(gene2bait_ref), "Gene2bait target:", len(gene2bait_target))
    print("Bait ref:", len(interaction_ref.keys()), "Bait target:", len(interaction_target.keys()))

    for gene_ref in ortho.keys():
        nb_contact = length_ref = 0
        # Are orthologous genes in baited fragments in both sides ?
        if gene_ref in gene2bait_ref.keys():
            if ortho[gene_ref] in gene2bait_target.keys():
                for bait_ref in gene2bait_ref[gene_ref]:
                    # Does reference bait in PC-HiC data ?
                    if bait_ref in interaction_ref.keys():
                        CAGE_stats = [[0]*5]
                        ENCODE_stats = [[0]*5]
                        RoadMap_stats = [[0]*5]
                        GRO_seq_stats = [[0]*5]

                        # Adding contact information
                        for PIR_ref in interaction_ref[bait_ref]:
                            contact_ref = PIR_ref[0]
                            length_ref += PIR_ref[1]
                            nb_contact += 1

                            # Does this contact contain any enhancer ?
                            CAGE_stats.append(conserv_enh("CAGE", contact_ref, gene_ref))
                            if ref_sp == "human":
                                ENCODE_stats.append(conserv_enh("ENCODE", contact_ref, gene_ref))
                                RoadMap_stats.append(conserv_enh("RoadMap", contact_ref, gene_ref))
                                GRO_seq_stats.append(conserv_enh("GRO_seq", contact_ref, gene_ref))

                        CAGE_stats = [sum(x) for x in zip(*CAGE_stats)]

                        if ref_sp != "human":
                            summary_conserv[gene_ref] = list(
                                itertools.chain([nb_contact], [length_ref], CAGE_stats)) #, ENCODE_stats))
                        else:
                            ENCODE_stats = [sum(x) for x in zip(*ENCODE_stats)]
                            RoadMap_stats = [sum(x) for x in zip(*RoadMap_stats)]
                            GRO_seq_stats = [sum(x) for x in zip(*GRO_seq_stats)]
                            summary_conserv[gene_ref] = list(itertools.chain([nb_contact], [length_ref], CAGE_stats, ENCODE_stats, RoadMap_stats, GRO_seq_stats))


    ################################################ Writting output ###################################################
    output_file = "/conservation/" + ref_sp + data_ref + "_merged_to_" + target_sp + data_target + "_summary_enhancers_by_gene.txt"
    output_summary = open(path_result + output_file, 'w')
    if os.stat(path_result + output_file).st_size == 0:
        if ref_sp == "human":
            output_summary.write(
                "gene\tnb_contact\tcontacted_length\t"
                "nb_CAGE\tnb_CAGE_conserv\tnb_CAGE_overlap\tnb_CAGE_synt\tnb_CAGE_contact\t"
                "nb_ENCODE\tnb_ENCODE_conserv\tnb_ENCODE_overlap\tnb_ENCODE_synt\tnb_ENCODE_contact\t"
                "nb_RoadMap\tnb_RoadMap_conserv\tnb_RoadMap_overlap\tnb_RoadMap_synt\tnb_RoadMap_contact\t"
                "nb_GRO_seq\tnb_GRO_seq_conserv\tnb_GRO_seq_overlap\tnb_GRO_seq_synt\tnb_GRO_seq_contact\n")
        else:
            output_summary.write("gene\tnb_contact\tcontacted_length\t"
                                 "nb_CAGE\tnb_CAGE_conserv\tnb_CAGE_overlap\tnb_CAGE_synt\tnb_CAGE_contact\n")
                                 #"nb_ENCODE\tnb_ENCODE_conserv\tnb_ENCODE_overlap\tnb_ENCODE_synt\tnb_ENCODE_contact\n")

    for gene in summary_conserv.keys():
        output_summary.write(gene + '\t' + str("\t".join(str(x) for x in summary_conserv[gene])) + '\n')

    output_summary.close()


data1 = ["", "_simul"]  # or "" for observed data
data2 = ["", "_simul"]

for ref in data1:
    for target in data2:
        print("Origin sp:", ref_sp, "; data:", ref, "; to data :", target)
        running_all(ref, target)
