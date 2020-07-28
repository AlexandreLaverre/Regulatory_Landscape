
import os
import numpy as np

ref_sp = "mouse"

path = "/home/laverre/Data/Regulatory_landscape/result"
path_evol = path + "/Supplementary_dataset6_regulatory_landscape_evolution/" + ref_sp + "/"
path_contact = path + "/Supplementary_dataset4_genes_enhancers_contacts/" + ref_sp + "/"


############################################# Enhancers contacts ###################################################
def gene_enh_contact(enh_name, data):
    enh_infos = {}
    enh_stats = {}
    with open(path_contact + "gene_" + enh_name + "_enhancers_" + data + "_interactions.txt") as f1:
        first_line = f1.readline().strip("\n")
        first_line = first_line.split("\t")
        sample_name = first_line[6:]

        for i in f1.readlines():
            i = i.strip("\n")
            i = i.split("\t")

            bait = i[0]
            gene = i[1]
            enh = i[3]
            dist = float(i[4])
            score = float(i[5])
            sample = [float(x) for x in i[6:len(i)]]
            info = [[bait], [gene]]

            if enh not in enh_infos.keys():
                enh_infos[enh] = info
                enh_stats[enh] = [[dist], [score], sample]
            else:
                for index, e in enumerate(enh_infos[enh]):
                    enh_infos[enh][index] = list(set(enh_infos[enh][index] + info[index]))

                enh_stats[enh][0].append(dist)
                enh_stats[enh][1].append(score)
                enh_stats[enh][2] = list(np.nansum((enh_stats[enh][2], sample), axis=0))

    output_file = path_evol + "/enhancers_conservation/" + enh_name + "/" + data + "_stats.txt"
    output = open(output_file, 'w')
    if os.stat(output_file).st_size == 0:
       output.write("enh\tnb_bait\tnb_gene\tmed_dist\tmed_score\tnb_sample\t" + "\t".join(sample_name) + "\n")

    for enhancer in enh_infos.keys():
        nb_sample = str(len([x for x in enh_stats[enhancer][2] if x > 0]))
        enh_stats[enhancer][2] = [str(1) if x > 0 else str(0) for x in enh_stats[enhancer][2]]

        output.write(enhancer + "\t" + str(len(enh_infos[enhancer][0])) + "\t" + str(len(enh_infos[enhancer][1])) + "\t" +
                     str(np.median(enh_stats[enhancer][0])) + "\t" + str(np.median(enh_stats[enhancer][1])) + "\t" +
                     nb_sample + "\t" + "\t".join(enh_stats[enhancer][2]) + "\n")

    output.close()


datas = ["original", "simulated"]
enh_datas = ["CAGE", "ENCODE"]
if ref_sp == "human":
   enh_datas.extend(["GRO_seq", "RoadMap"])

for dat in datas:
    for enh_data in enh_datas:
        gene_enh_contact(enh_data, dat)
        print(enh_data, "in", dat, "data done !")

