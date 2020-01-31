#!/usr/bin/env python3
# coding=utf-8

import os
import numpy as np

# Conservation mouse interaction in human:
ref_sp = "human"
target_sp = "mouse"
data_ref = ""  # or "" for observed data
data_target = ""

path_data = "/home/laverre/Documents/Regulatory_Landscape/data/"
path_result = "/home/laverre/Documents/Regulatory_Landscape/result/"
seuil = 0.1

print("Origin sp:", ref_sp, "; data:", data_ref, "; to data :", data_target)


# Chromatin contact dictionary in both species
def interaction_dict(sp, data):
    interaction = {}
    stats = {}
    if data == "_simul":
        input = "/Simulations/simulations_" + sp + "_10Mb_bin5kb_fragoverbin_chr_merged.txt"
    else:
        input = "/all_interactions/all_interactions_chr_merged.txt_cell_names"

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
                stats[bait + '-' + PIR] = (dist, nb_type, median_strength)

                if 25000 < abs(dist) <= 10000000:
                    if bait in interaction.keys():
                        interaction[bait].append(PIR)
                    else:
                        interaction[bait] = [PIR]

    return interaction


interaction_ref = interaction_dict(ref_sp, data_ref)
interaction_target = interaction_dict(target_sp, data_target)
print("Interaction dictionary done !")


# Genes to Bait
def dict_gene2bait(sp, data):
    gene2bait = {}
    with open(path_result + "conservation/bait_composition_" + sp + data + "_merged.txt") as f3:
        for i in f3.readlines()[1:]:
            i = i.strip("\n")
            i = i.split("\t")
            bait = (str(i[0]) + ":" + str(i[1]) + ":" + str(i[2]))

            if i[len(i)-1] != '':
                genes = i[len(i)-1].split(',')
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
with open(
        path_result + "conservation/alignments/" + ref_sp + "/contacted_seq/AlignmentStatistics_Excluding_all_Exons_" + ref_sp + "2" + target_sp + "_merged_bait.txt") as f1:
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
        overlap_target[frag_align] = i[4].split(',')  # Add all overlapping fragment

print("Overlap aligned with target sp frag : done !")

# Duplication score
frag_dupli = {}
with open(path_result + "BLAT_duplication/" + ref_sp + "_merged_fragments_duplication_stats.txt") as f1:
    for i in f1.readlines()[1:]:
        i = i.strip("\n")
        i = i.split("\t")
        frag = i[0]
        frag_dupli[frag] = int(i[3]) - 1

print("Score dupli : done ! ")


print("Running conservation of interactions... ")
conserv_contact = {}
for gene_ref in ortho.keys():
    # Are orthologous gene in baited fragments in both side
    if gene_ref in gene2bait_ref.keys():
        if ortho[gene_ref] in gene2bait_target.keys():
            for bait_ref in gene2bait_ref[gene_ref]:
                for bait_target in gene2bait_target[ortho[gene_ref]]:

                    # Are baited fragment in PC-HiC data in both side
                    if bait_ref in interaction_ref.keys():
                        if bait_target in interaction_target.keys():

                            # Is contacted fragment conserved in target sp
                            for contact_ref in interaction_ref[bait_ref]:
                                if contact_ref in frag_conserv.keys():
                                    conserv_frag = frag_conserv[contact_ref][0]
                                    overlap_score = frag_conserv[contact_ref][1]

                                    # Is this fragment overlap with contacted fragment in target sp
                                    if conserv_frag in overlap_target.keys():
                                        for homolog_frag in overlap_target[conserv_frag]:

                                            # Is this homolog fragment in contact with ortholog gene
                                            if homolog_frag in interaction_target[bait_target]:
                                                conserv_contact[gene_ref + '\t' + contact_ref] \
                                                    = (gene_target, homolog_frag, overlap_score)

                                            elif gene_ref + '\t' + contact_ref not in conserv_contact.keys():
                                                conserv_contact[gene_ref + '\t' + contact_ref] \
                                                    = ("not_good_gene", homolog_frag, overlap_score)

                                    elif gene_ref + '\t' + contact_ref not in conserv_contact.keys():
                                        conserv_contact[gene_ref + '\t' + contact_ref] \
                                            = ("not_in_contact", conserv_frag, overlap_score)

                                elif gene_ref + '\t' + contact_ref not in conserv_contact.keys():
                                    conserv_contact[gene_ref + '\t' + contact_ref] = ("not_conserved", "NA", "NA")


output_file = "/conservation/interaction_conservation/" + ref_sp + data_ref + "_merged_to_" + target_sp + data_target + "_by_gene.txt"

output = open(path_result + output_file, 'w')
if os.stat(path_result + output_file).st_size == 0:
    output.write("gene\tcontact\tgene_target\thomologous_contact\toverlap\tduplication\tconserv\n")

for inter in conserv_contact.keys():
    gene_ref = inter.split('\t')[0]
    contact_ref = inter.split('\t')[1]

    if contact_ref not in frag_dupli.keys():
        frag_dupli[contact_ref] = 0

    if contact_ref not in frag_conserv.keys():
        frag_conserv[contact_ref] = ("NA", 0)

    output.write(gene_ref + '\t' + contact_ref + '\t' + str(conserv_contact[inter][0]) + '\t'
                 + str(conserv_contact[inter][1]) + '\t' + str(conserv_contact[inter][2]) + '\t'
                 + str(frag_dupli[contact_ref]) + '\t' + str(frag_conserv[contact_ref][0]) + '\n')

output.close()
