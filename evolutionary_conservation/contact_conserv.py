#!/usr/bin/env python3
# coding=utf-8

import os

ref_sp = "mouse"
target_sp = "human" if ref_sp == "mouse" else "mouse"

path = "/home/laverre/Data/Regulatory_landscape/result"
path_annot = path + "/Supplementary_dataset3_annotations/" + ref_sp + "/"
path_evol = path + "/Supplementary_dataset6_regulatory_landscape_evolution/" + ref_sp + "/"
path_contact = path + "/Supplementary_dataset4_genes_enhancers_contacts/"


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


###################################### Enhancers statistics & Alignments ##############################################
def enh_info(enh_name):
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

    return align, stats


############################## Contact map in target species with lifted enhancers ####################################
def gene_enh_target(data, enh):
    contact = {}
    stats = {}
    enh_name = "lifted_" + enh
    with open(path_contact + target_sp + "/gene_" + enh_name + "_enhancers_" + data + "_interactions.txt") as f1:
        first_line = f1.readline().strip("\n")
        first_line = first_line.split("\t")
        sample_name = first_line[7:]

        for i in f1.readlines():
            i = i.strip("\n")
            i = i.split("\t")

            gene = i[1]
            enh = i[3]
            samples = [str(x) if x != 'nan' else str(0) for x in i[7:len(i)]]
            nb_sample = str(len([x for x in i[7:len(i)] if x != 'nan']))
            infos = [[i[4]], [i[5]], [nb_sample], samples]
            stats[gene + '-' + enh] = [e for sublist in infos for e in sublist]  # dist, median score, samples

            if gene not in contact.keys():
                contact[gene] = [enh]
            else:
                contact[gene].append(enh)

    return contact, stats, sample_name


################################### Conservation of contacts between species #########################################
def conserv_contact(data, enh_name):
    output_file = path_evol + "contact_conservation/" + enh_name + "/" + ref_sp + "2" + target_sp + "_" + data + ".txt2"
    output = open(output_file, 'w')

    conserv_enh, stats_enh = enh_info(enh_name)

    target_contact, target_stats, target_sample_name = gene_enh_target(data, enh_name)

    with open(path_contact + ref_sp + "/gene_" + enh_name + "_enhancers_" + data + "_interactions.txt") as f1:
        first_line = f1.readline().strip("\n")
        first_line = first_line.split("\t")
        sample_name = first_line[7:]

        if os.stat(output_file).st_size == 0:
            output.write("origin_gene\torigin_enh\torigin_dist\torigin_med_score\tnb_sample\t" + "\t".join(sample_name) + "\t"
                         "length_enh\trepeat_part\tGC_rate\tBLAT_match\talign_score\t"
                         "target_gene\tlifted_enh\ttarget_dist\ttarget_med_score\ttarget_nb_sample\t" + "\t".join(target_sample_name)+"\n")

        for i in f1.readlines():
            i = i.strip("\n")
            i = i.split("\t")

            gene = i[1]
            enh = i[3]
            samples = [str(x) if x != 'nan' else str(0) for x in i[7:len(i)]]
            nb_sample = str(len([x for x in i[7:len(i)] if x != 'nan']))

            # Do orthologous genes in PC-HiC data in both sides ?
            if gene in ortho.keys() and ortho[gene] in target_contact.keys():

                # Do enhancer conserved in target sp ?
                if enh in conserv_enh.keys():
                    lifted_enh = conserv_enh[enh][0]
                    align_score = str(conserv_enh[enh][1])

                    # Do conserved enhancer in contact with orthologous gene in target sp ?
                    if lifted_enh in target_contact[ortho[gene]]:
                        output.write(gene + "\t" + enh + "\t" + i[4] + "\t" + i[5] + "\t" + nb_sample + '\t')
                        output.write("\t".join(samples) + "\t")
                        output.write("\t".join(stats_enh[enh]) + "\t" + align_score + "\t")

                        output.write(ortho[gene] + "\t" + lifted_enh + "\t" +
                                     "\t".join(target_stats[ortho[gene]+'-'+lifted_enh]) + "\n")

    output.close()


datas = ["simulated"]
enhancers = ["CAGE", "ENCODE"]
if ref_sp == "human":
     enhancers.extend(["GRO_seq", "RoadMap"])


for dat in datas:
    print("Running", ref_sp, "to", target_sp, "in", dat, "contacts :")
    for enhancer in enhancers:
        conserv_contact(dat, enhancer)
        print(enhancer, "done !")



