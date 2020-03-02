#!/usr/bin/env python3
# coding=utf-8

import os

ref_sp = "mouse"

path = "/home/laverre/Data/Regulatory_landscape/result"
path_annot = path + "/Supplementary_dataset3_annotations/" + ref_sp + "/"
path_evol = path + "/Supplementary_dataset6_regulatory_landscape_evolution/" + ref_sp + "/"
path_contact = path + "/Supplementary_dataset4_genes_enhancers_contacts/" + ref_sp + "/"


############################################### Enhancers statistics ##################################################
def stats_enh(enh_name):
    stats = {}
    with open(path_annot + enh_name + "_BLAT_summary_0.8.txt") as f1:
        for i in f1.readlines()[1:]:
            i = i.strip("\n")
            i = i.split("\t")

            enh = i[0]
            repeat_part = int(i[2]) / int(i[1])  # repeat = nb_N / length
            try:
                GC_rate = int(i[3]) / (int(i[1]) - int(i[2]))  # GC = nb_GC / (length - nb_N)
            except ZeroDivisionError:
                GC_rate = "NA"

            # enh = (length, unrepeat, GC, nb_match)
            stats[enh] = [str(i[1]), str(repeat_part), str(GC_rate), str(i[4])]

    return stats


def synt_conserv(target_sp, data):
    ########################################### Orthologous gene one2one ###############################################
    ortho = {}
    origin_gene_coord = {}
    target_gene_coord = {}
    with open(path_evol + "gene_orthology/" + ref_sp + "2" + target_sp + "_orthologue.txt") as f3:
        for i in f3.readlines()[1:]:
            i = i.strip("\n")
            i = i.split("\t")
            if i[9] == "ortholog_one2one":
                origin_gene_coord[str(i[0])] = ["chr" + str(i[1]), i[2], i[3]]  # ID = (chr, start, end)
                target_gene_coord[str(i[5])] = ["chr" + str(i[6]), i[7], i[8]]

                ortho[str(i[0])] = str(i[5])

    ############################################ Enhancers alignments ###############################################
    def align_enh(enh_name):
        align = {}
        with open(path_evol + "enhancers_conservation/" + enh_name + "/AlignmentStatistics_Excluding_all_Exons_" +
                  ref_sp + "2" + target_sp + "_" + enh_name + ".txt") as f1:
            for i in f1.readlines()[1:]:
                i = i.strip("\n")
                i = i.split("\t")

                enh_origin = i[0].strip(":+")
                enh_target = i[1].strip(":+")
                enh_target = enh_target.strip(":-")

                try:
                    align_score = int(i[6]) / int(i[9])  # FilteredUngappedLength / FilteredAlignmentLength
                except ZeroDivisionError:
                    align_score = 0

                align[enh_origin] = (enh_target, align_score)

        return align

    ############################################ Synteny conserv ###############################################
    def gene_enh_contact(enh_name):
        output_file = path_evol + "synteny_conservation/" + enh_name + "/" + ref_sp + "2" + target_sp \
                      + "_" + enh_name + "_" + data + "_synteny.txt"
        output = open(output_file, 'w')
        if os.stat(output_file).st_size == 0:
            output.write("origin_gene\torigin_gene_coord\torigin_enh\torigin_dist\t"
                         "target_gene\ttarget_gene_coord\ttarget_enh\ttarget_dist\t"
                         "length_enh\trepeat_part\tGC_rate\tBLAT_match\n")

        dic_aligned = align_enh(enh_name)
        dic_stats = stats_enh(enh_name)

        with open(path_contact + "gene_" + enh_name + "_enhancers_" + data + "_interactions.txt") as f1:
            for i in f1.readlines()[1:]:
                i = i.strip("\n")
                i = i.split("\t")

                origin_gene = i[1]
                origin_enh = i[3]
                origin_dist = i[4]

                if origin_gene in ortho.keys():
                    if origin_enh in dic_aligned.keys():

                        target_gene = ortho[origin_gene]
                        target_gene_chr = target_gene_coord[target_gene][0]
                        target_TSS = target_gene_coord[target_gene][1]

                        target_enh = dic_aligned[origin_enh][0]
                        chr_target_enh = target_enh.split(':')[0]
                        start_target_enh = target_enh.split(':')[1]
                        end_target_enh = target_enh.split(':')[2]
                        mid_target_enh = (int(start_target_enh)+int(end_target_enh)) / 2

                        if target_gene_chr == chr_target_enh:
                            target_dist = abs(int(target_TSS) - mid_target_enh)
                        else:
                            target_dist = "trans"

                        output.write(origin_gene + "\t" + ':'.join(origin_gene_coord[origin_gene]) + "\t" +
                                     origin_enh + "\t" + origin_dist + "\t" +
                                     target_gene + "\t" + ':'.join(target_gene_coord[target_gene]) + "\t" +
                                     target_enh + "\t" + str(target_dist) + "\t" + '\t'.join(dic_stats[origin_enh]) + "\n")

        output.close()

    print("Running ENCODE...")
    gene_enh_contact("ENCODE")
    print("Running CAGE...")
    gene_enh_contact("CAGE")
    if ref_sp == "human":
        print("Running RoadMap...")
        gene_enh_contact("RoadMap")
        print("Running GRO_seq...")
        gene_enh_contact("GRO_seq")


datas = ["original", "simulated"]
species = ["human", "cow", "opossum", "elephant", "rabbit", "rat", "macaque", "dog", "chicken"]

for dat in datas:
    for sp in species:
        print("Running ", ref_sp, "to", sp, "in", dat, "contacts :")
        synt_conserv(sp, dat)






