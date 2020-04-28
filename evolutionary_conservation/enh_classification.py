#!/usr/bin/env python3
# coding=utf-8

ref_sp = "mouse"
path = "/home/laverre//Documents/Regulatory_Landscape/data/"


############################################### Enhancers informations #################################################
def overlap_enh(data, file):
    overlap = {}
    with open(path + ref_sp + "/overlap/" + data + file) as f1:
        for i in f1.readlines()[1:]:
            i = i.strip("\n")
            i = i.split("\t")

            enh = i[0].strip(":+")
            enh = enh.strip(":-")
            nb_overlap_bp = int(i[6]) if file == "_overlap_all_exons.txt" or file == "_overlap_all_genes.txt" else int(i[4])

            overlap[enh] = nb_overlap_bp

    return overlap


enhancers = ["CAGE", "ENCODE", "RoadMap", "GRO_seq"] if ref_sp == "human" else ["CAGE"]

for enh_name in enhancers:
    ov_exon = overlap_enh(enh_name, "_overlap_all_exons.txt")
    ov_gene = overlap_enh(enh_name, "_overlap_all_genes.txt")
    ov_genic_50k = overlap_enh(enh_name, "_overlap_all_genes_50kb.txt")
    ov_genic_100k = overlap_enh(enh_name, "_overlap_all_genes_100kb.txt")
    ov_genic_500k = overlap_enh(enh_name, "_overlap_all_genes_500kb.txt")
    ov_exonic_50k = overlap_enh(enh_name, "_overlap_all_exons_50kb.txt")
    ov_exonic_100k = overlap_enh(enh_name, "_overlap_all_exons_100kb.txt")
    ov_exonic_500k = overlap_enh(enh_name, "_overlap_all_exons_500kb.txt")

    output_file = path + ref_sp + "/" + enh_name + "_classification.txt"
    output = open(output_file, 'w')
    output.write("enh\tstart\tend\ttype\t50kb_genic_bp\t100kb_genic_bp\t500kb_genic_bp\t"
                 "50kb_exonic_bp\t100kb_exonic_bp\t500kb_exonic_bp\n")

    for enhancer in ov_exon.keys():
        start = str(enhancer.split(':')[1])
        end = str(enhancer.split(':')[2])

        # Classification
        if ov_exon[enhancer] > 0:
            type = "exonic"
        elif ov_gene[enhancer] > 0:
            type = "genic"
        else:
            type = "intergenic"

        output.write(enhancer + '\t' + start + '\t' + end + '\t' + type + '\t')
        output.write(str(ov_genic_50k[enhancer]) + '\t' + str(ov_genic_100k[enhancer]) + '\t' + str(ov_genic_500k[enhancer]) + '\t')
        output.write(str(ov_exonic_50k[enhancer]) + '\t' + str(ov_exonic_100k[enhancer]) + '\t' + str(ov_exonic_500k[enhancer]) + '\n')

    output.close()

