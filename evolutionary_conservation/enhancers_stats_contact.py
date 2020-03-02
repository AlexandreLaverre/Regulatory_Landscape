
import os
import numpy as np

ref_sp = "human"

path = "/home/laverre/Data/Regulatory_landscape/result"
path_evol = path + "/Supplementary_dataset6_regulatory_landscape_evolution/" + ref_sp + "/"
path_contact = path + "/Supplementary_dataset4_genes_enhancers_contacts/" + ref_sp + "/"


############################################# Enhancers contacts ###################################################
def gene_enh_contact(enh_name, data):
    enh_infos = {}
    enh_stats = {}
    with open(path_contact + "gene_" + enh_name + "_enhancers_" + data + "_interactions.txt") as f1:
        for i in f1.readlines()[1:]:
            i = i.strip("\n")
            i = i.split("\t")

            bait = i[0]
            gene = i[1]
            enh = i[3]
            dist = float(i[4])
            score = float(i[5])
            sample = i[7].split(',')
            info = [[bait], [gene], [i for i in sample]]

            if enh not in enh_infos.keys():
                enh_infos[enh] = [[bait], [gene], [i for i in sample]]
                enh_stats[enh] = [[dist], [score]]
            else:
                for index, e in enumerate(enh_infos[enh]):
                    enh_infos[enh][index] = list(set(enh_infos[enh][index] + info[index]))

                enh_stats[enh][0].append(dist)
                enh_stats[enh][1].append(score)

    output_file = path_evol + enh_name +    "/enhancers_conservation/" + enh_name + "_" + data + "_stats.txt"
    output = open(output_file, 'w')
    if os.stat(output_file).st_size == 0:
        output.write("enh\tnb_bait\tnb_gene\tsamples\tmed_dist\tmed_score\n")

    for enhancer in enh_infos.keys():
        output.write(enhancer + "\t" + str(len(enh_infos[enhancer][0])) + "\t" + str(len(enh_infos[enhancer][1])) + "\t" +
                     ",".join(enh_infos[enhancer][2]) + "\t" + str(np.median(enh_stats[enhancer][0])) + "\t" +
                     str(np.median(enh_stats[enhancer][1])) + "\n")

    output.close()


datas = ["original", "simulated"]
enh_datas = ["CAGE", "ENCODE"]
if ref_sp == "human":
    enh_datas.extend(["GRO_seq", "RoadMap"])

for dat in datas:
    for enh_data in enh_datas:
        gene_enh_contact(enh_data, dat)
        print(enh_data, "in", dat, "data done !")

