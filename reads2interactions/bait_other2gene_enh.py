#!/usr/bin/env python3
# coding=utf-8

import os
import numpy as np

ref_sp = "mouse"
target_sp = "human" if ref_sp == "mouse" else "mouse"

path = "/home/laverre/Data/Regulatory_landscape/result"
path_annot = path + "/Supplementary_dataset3_annotations/"
path_result = path + "/Supplementary_dataset4_genes_enhancers_contacts/" + ref_sp + "/"
path_gene = path + "/Supplementary_dataset6_regulatory_landscape_evolution/" + ref_sp + "/"


#################################################### Genes TSS ########################################################
gene_TSS = {}
with open(path_annot + "genes/" + ref_sp + "_genes_Ensembl94.txt") as f3:
    for i in f3.readlines()[1:]:
        i = i.strip("\n")
        i = i.split("\t")

        gene_TSS[str(i[5])] = i[4]  # ID = (chr, start, end)

#################################################### Bait to Genes ####################################################
bait2gene = {}
with open(path_annot + ref_sp + "/bait_overlap_TSS_1kb.txt") as f1:
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


########################################### Contacted Fragment to Enhancers ###########################################
def fragment2enhancer(enh_data):
    frag2enh = {}
    with open(path_annot + ref_sp + "/" + enh_data) as f1:
        for i in f1.readlines()[1:]:
            i = i.strip("\n")
            i = i.split("\t")
            frag = (str(i[1]) + ":" + str(i[2]) + ":" + str(i[3]))

            if i[4] != "NA":
                enhancers = i[4].split(',')
                for enh in enhancers:
                    if frag not in frag2enh.keys():
                        frag2enh[frag] = [enh]
                    else:
                        frag2enh[frag].append(enh)

    return frag2enh


ENCODE = fragment2enhancer("frag_overlap_ENCODE.txt")
lifted_ENCODE = fragment2enhancer("frag_overlap_lifted_"+target_sp+"_ENCODE.txt")
CAGE = fragment2enhancer("frag_overlap_CAGE.txt")
lifted_CAGE = fragment2enhancer("frag_overlap_lifted_"+target_sp+"_CAGE.txt")

if ref_sp == "human":
    RoadMap = fragment2enhancer("frag_overlap_RoadMap.txt")
    GRO_seq = fragment2enhancer("frag_overlap_GRO_seq.txt")

if ref_sp == "mouse":
    lifted_RoadMap = fragment2enhancer("frag_overlap_lifted_"+target_sp+"_RoadMap.txt")
    lifted_GRO_seq = fragment2enhancer("frag_overlap_lifted_"+target_sp+"_GRO_seq.txt")


####################################### PC-HIC Contacts to Regulatory Landscape #######################################
def gene_enh_contact(data, data_name, enh_data, enh_name):
    output_file = path_result + "gene_" + enh_name + "_enhancers_" + data_name + "_interactions.txt"
    output = open(output_file, 'w')
    if os.stat(output_file).st_size == 0:
        output.write("bait\tgene\tcontact\tenhancer\tdist\tmedian_score\tnb_sample\tsamples\n")

    with open(data) as f1:
        first_line = f1.readline().strip("\n")
        first_line = first_line.split("\t")
        sample_name = first_line[8:]

        for i in f1.readlines():
            i = i.strip("\n")
            i = i.split("\t")

            bait = (str(i[0]) + ":" + str(i[1]) + ":" + str(i[2]))
            contacted = (str(i[3]) + ":" + str(i[4]) + ":" + str(i[5]))
            contact_score = [float(x) for x in i[8:] if x != "NA"]
            median_score = str(np.median(contact_score))
            contact_sample = [sample_name[idx] for idx, x in enumerate(i[8:]) if x != "NA"]
            nb_sample = str(len(contact_score))

            if i[0] == i[3]:    # cis-interaction
                dist = float(i[7])
                if 25000 < dist < 10000000:
                    if i[6] == "unbaited":      # bait-other interaction

                        if bait in bait2gene.keys():    # Is bait contain gene ?
                            for gene in bait2gene[bait]:
                                if contacted in enh_data.keys():    # Is contacted fragment contain enhancer ?
                                    for enh in enh_data[contacted]:

                                        enh_coord = enh.split(":")
                                        mid_enh = (int(enh_coord[2]) + int(enh_coord[1])) / 2
                                        dist_gene_enh = abs(float(gene_TSS[gene])-mid_enh)

                                        output.write(bait + "\t" + gene + "\t" + contacted + "\t" + enh + "\t")
                                        output.write(str(dist_gene_enh) + "\t" + median_score + "\t" + nb_sample + "\t")
                                        output.write(str(','.join(contact_sample)) + "\n")

    output.close()


data_obs = path + "/Supplementary_dataset1_original_interactions/" + ref_sp + "/all_interactions.txt"
data_simul = path + "/Supplementary_dataset2_simulated_interactions/" + ref_sp + "/simulated_all_interactions.txt"
datas = [data_obs, data_simul]

print("Running script on ", ref_sp)

for dat in datas:
    name = "original" if dat == data_obs else "simulated"

    #gene_enh_contact(dat, name, ENCODE, "ENCODE")
    print(name, "ENCODE done !")
    gene_enh_contact(dat, name, lifted_ENCODE, "lifted_ENCODE")
    print(name, "lifted ENCODE done !")
    gene_enh_contact(dat, name, CAGE, "CAGE")
    print(name, "CAGE done !")
    gene_enh_contact(dat, name, lifted_CAGE, "lifted_CAGE")
    print(name, "lifted CAGE done !")

    if ref_sp == "human":
        gene_enh_contact(dat, name, RoadMap, "RoadMap")
        print(name, "RoadMap done !")
        gene_enh_contact(dat, name, GRO_seq, "GRO_seq")
        print(name, "GRO_seq done !")

    if ref_sp == "mouse":
        gene_enh_contact(dat, name, lifted_RoadMap, "lifted_RoadMap")
        print(name, "lifted RoadMap done !")
        gene_enh_contact(dat, name, lifted_GRO_seq, "lifted_GRO_seq")
        print(name, "lifted GRO_seq done !")




