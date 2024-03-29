#!/usr/bin/env python3
# coding=utf-8

import os

ref_sp = "mouse"
target_sp = "human" if ref_sp == "mouse" else "mouse"

path = "/home/laverre/Manuscript"
path_annot = path + "/SupplementaryDataset4/" + ref_sp + "/"
path_evol = path + "/SupplementaryDataset7/" + ref_sp + "/"
path_contact = path + "/SupplementaryDataset4/"


########################################### Orthologous gene one2one ###############################################
ortho = {}

with open(path_evol + "gene_orthology/" + ref_sp + "2" + target_sp + "_orthologue.txt") as f3:
    for i in f3.readlines()[1:]:
        i = i.strip("\n")
        i = i.split("\t")
        if i[9] == "ortholog_one2one":
            ortho[str(i[0])] = str(i[5])

target_gene_TSS = {}
with open(path + "/SupplementaryDataset3/genes/" + target_sp + "_genes_Ensembl94.txt") as f3:
    for i in f3.readlines()[1:]:
        i = i.strip("\n")
        i = i.split("\t")
        TSS = str(i[3]) if i[5] == "1" else str(i[4])
        chr = "chr" + str(i[2]) if str(i[2]) != "MT" else "chrM"

        target_gene_TSS[str(i[0])] = [chr, TSS]  # ID = [chr, TSS]

###################################### Enhancers statistics & Alignments ##############################################
def enh_info(enh_name, data):
    stats = {}
    with open(path_annot + enh_name + "/statistics_contacted_enhancers_" + data + ".txt") as f1:
        for i in f1.readlines()[1:]:
            i = i.strip("\n")
            i = i.split("\t")

            enh = i[0]+':'+i[1]+':'+i[2]
            length = int(i[3])
            repeat_part = int(i[13]) / length  # repeat = nb_N / length
            try:
                GC_rate = int(i[14]) / (length - int(i[13]))  # GC = nb_GC / (length - nb_N)
            except ZeroDivisionError:
                GC_rate = "NA"

            # enh = (length, unrepeat, GC, nb_match)
            stats[enh] = [str(length), str(repeat_part), str(GC_rate), str(i[11])]

    align = {}
    with open(path_evol + "sequence_conservation/enhancers/" + enh_name + "/AlignmentStatistics_Excluding_Exons_" +
              ref_sp + "2" + target_sp + ".txt") as f1:
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
def gene_enh_target(data_target, enh):
    contact = {}
    stats = {}
    with open(path_contact + target_sp + "/" + enh + "/gene_lifted_enhancer_contacts_" + data_target + "_interactions.txt") as f1:
        first_line = f1.readline().strip("\n")
        first_line = first_line.split("\t")
        sample_name = first_line[7:]

        for i in f1.readlines():
            i = i.strip("\n")
            i = i.split("\t")

            gene = i[0]
            enh = i[1]
            samples = [str(x) if x != 'nan' else str(0) for x in i[7:len(i)]]
            nb_sample = str(len([x for x in i[7:len(i)] if x != 'nan']))
            infos = [[i[5]], [nb_sample], samples]      # dist="i[4]
            stats[gene + '-' + enh] = [e for sublist in infos for e in sublist]  # median score, samples

            if gene not in contact.keys():
                contact[gene] = [enh]
            else:
                contact[gene].append(enh)

    return contact, stats, sample_name


################################### Conservation of contacts between species #########################################
def conserv_contact(data, data_target, enh_name):
    output_file = path_evol + "contact_conservation/" + enh_name + "/" + ref_sp + "_" + data + "2" + target_sp + "_" + data_target + ".txt_test"
    output = open(output_file, 'w')

    conserv_enh, stats_enh = enh_info(enh_name, data)

    target_contact, target_stats, target_sample_name = gene_enh_target(data_target, enh_name)

    with open(path_contact + ref_sp + "/" + enh_name + "/gene_enhancer_contacts_" + data + "_interactions.txt") as f1:
        first_line = f1.readline().strip("\n")
        first_line = first_line.split("\t")
        sample_name = first_line[7:]

        if os.stat(output_file).st_size == 0:
            output.write("origin_gene\torigin_enh\torigin_dist\torigin_med_score\tnb_sample\t" + "\t".join(sample_name) + "\t"
                         "length_enh\trepeat_part\tGC_rate\tBLAT_match\ttarget_gene\ttarget_data\t"
                         "lifted_enh\talign_score\ttarget_dist\ttarget_med_score\ttarget_nb_sample\t" + "\t".join(target_sample_name)+"\n")

        for i in f1.readlines():
            i = i.strip("\n")
            i = i.split("\t")

            gene = i[0]
            enh = i[1]
            samples = [str(x) if x != 'nan' else str(0) for x in i[7:len(i)]]
            nb_sample = str(len([x for x in i[7:len(i)] if x != 'nan']))

            output.write(gene + "\t" + enh + "\t" + i[4] + "\t" + i[5] + "\t" + nb_sample + '\t')
            output.write("\t".join(samples) + "\t" + "\t".join(stats_enh[enh]) + "\t")

            target_gene = lifted_enh = target_dist = target_gene_in_data = "NA"
            align_score = 0
            conserved_stats = "\tNA" * 2 + "\t0" * len(target_sample_name)

            # Do enhancer conserved in target sp ?
            if enh in conserv_enh.keys():
                lifted_enh = conserv_enh[enh][0]
                align_score = conserv_enh[enh][1]

            # Do genes are orthologous ?
            if gene in ortho.keys():
                target_gene = ortho[gene]

                # Calculate distance in target
                if lifted_enh != "NA":
                    if target_gene_TSS[target_gene][0] == lifted_enh.split(":")[0]:   # same chromosome
                        mid_enh = (int(lifted_enh.split(":")[1])+int(lifted_enh.split(":")[2]))/2
                        target_dist = abs(int(target_gene_TSS[target_gene][1])-mid_enh)
                    else:
                        target_dist = "trans"

                # Do orthologous genes in PC-HiC data in both sides ?
                if target_gene in target_contact.keys():
                    target_gene_in_data = "TRUE"
                    # Do conserved enhancer in contact with orthologous gene in target sp ?
                    if lifted_enh in target_contact[target_gene]:
                        conserved_stats = "\t" + "\t".join(target_stats[target_gene+'-'+lifted_enh])

            output.write(target_gene + "\t" + target_gene_in_data + "\t" + lifted_enh + "\t" + str(align_score)
                         + "\t" + str(target_dist) + conserved_stats + "\n")

    output.close()


datas = ["original", "simulated"]
datas_target = ["original", "simulated"]
enhancers = ["FANTOM5", "ENCODE"]
if ref_sp == "human":
    enhancers.extend(["RoadmapEpigenomics", "FOCS_GRO_seq"])

for dat in datas:
    for dat_target in datas_target:
        print("Running", ref_sp, "in", dat, "to", target_sp, "in", dat_target, "contacts :")
        for enhancer in enhancers:
            conserv_contact(dat, dat_target, enhancer)
            print(enhancer, "done !")



