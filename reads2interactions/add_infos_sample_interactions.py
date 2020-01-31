#!/usr/bin/env python3
# coding=utf-8

import os
import glob

# Conservation mouse interaction in human:
sp = "mouse"
data = "_simul"  # or "_simul"
path_data = "/home/laverre/Documents/Regulatory_Landscape/data/" + sp + "/"

# path_data = "/beegfs/data/alaverre/Regulatory_landscape/result/simulations/"+sp+"_samples/"

if data == "_simul":
    data_type = "simulated"
    # all_interactions = path_data + "/Simulations/simulations_" + sp + "_10Mb_bin5kb_fragoverbin_chr.txt"
    all_interactions = path_data + "simulated_all_interactions.txt"
else:
    data_type = "observed"
    all_interactions = path_data + "/all_interactions/all_interactions_chr.txt"

if sp == "mouse":
    genome = "mm10"
else:
    genome = "hg38"

print("specie:", sp)
print("data:", data_type)

# PIR overlap bait
PIR_bait = {}
with open(
        "/beegfs/data/alaverre/Regulatory_landscape/data/CHICAGO_files/" + sp + "/" + genome + "/Digest_" + genome + "_HindIII_None.txt.baitmap") as f1:
    for i in f1.readlines():
        i = i.strip("\n")
        i = i.split("\t")
        frag = str(i[0]) + "\t" + str(i[1]) + "\t" + str(i[2])
        PIR_bait[frag] = "bait"

# All interactions dic --> nb_type info
all_dic = {}
with open(all_interactions, 'r') as f:
    for i in f.readlines()[1:]:
        i = i.strip("\n")
        i = i.split("\t")
        bait = (str(i[0]) + "\t" + str(i[1]) + "\t" + str(i[2]))
        PIR = (str(i[3]) + "\t" + str(i[4]) + "\t" + str(i[5]))

        if PIR in PIR_bait.keys():
            PIR_status = "baited"
        else:
            PIR_status = "unbaited"

        # Nb cell type
        if sp == "mouse":
            types = [i[6], min(i[7:11]), min(i[11], i[12]), i[13], min(i[14], i[15]), i[16], min(i[17:20])]
            nb_type = len([float(x) for x in types if x != "NA"])

        elif sp == "human":
            types = [min(i[6], i[22], i[30]), i[7], min(i[8:11]), i[11], i[12], i[13], i[14], i[15], i[16],
                     i[17], min(i[18], i[31]), i[19], i[20], i[21], min(i[23], i[24], i[25], i[29]),
                     min(i[26], i[27], i[28])]
            nb_type = len([float(x) for x in types if x != "NA"])

        all_dic[bait + ':' + PIR] = (nb_type, PIR_status)


# Sample interaction dic
def sample_dictionary(sample):
    sample_dic = {}
    with open(sample, 'r') as f:
        for i in f.readlines()[1:]:
            i = i.strip("\n")
            i = i.split("\t")
            bait = (str(i[0]) + "\t" + str(i[1]) + "\t" + str(i[2]))  # chr, start, end # +"chr"

            if data == "_simul":
                PIR = (str(i[3]) + "\t" + str(i[4]) + "\t" + str(i[5]), "NA", "NA")  # +"chr"
            else:
                PIR = (str(i[3]) + "\t" + str(i[4]) + "\t" + str(i[5]), str(i[6]),
                       str(i[7]))  # chr, start, end, N_reads, score # +"chr"

            if bait not in sample_dic.keys():
                sample_dic[bait] = [PIR]
            else:
                sample_dic[bait].append(PIR)

    # Writting output
    output_file = sample + "_infos"
    output = open(output_file, 'w')
    if os.stat(output_file).st_size == 0:
        output.write(
            "bait_chr\tbait_start\tbait_end\tchr\tstart\tend\tN_reads\tscore\tbaited_frag\tdist\tnb_contact\tnb_type\n")

    for bait, contacted in sample_dic.items():
        for i in contacted:
            if bait.split('\t')[0] == i[0].split('\t')[0]:
                midbait = ((int(bait.split('\t')[1]) + int(bait.split('\t')[2])) / 2)
                midcontact = ((int(i[0].split('\t')[1]) + int(i[0].split('\t')[2])) / 2)
                dist_obs = abs(midbait - midcontact)
            else:
                dist_obs = "NA"

            PIR = str(i[0])
            output.write(
                bait + "\t" + PIR + "\t" + i[1] + "\t" + i[2] + "\t" + str(str(all_dic[bait + ':' + PIR][1])) + "\t" +
                str(dist_obs) + "\t" + str(len(sample_dic[bait])) + "\t" + str(all_dic[bait + ':' + PIR][0]) + "\n")

    output.close()


if data == "_simul":
    samples = glob.glob(path_data + "*txt")
else:
    samples = glob.glob(path_data + "all_interactions/" + sp + "_samples/*ibed")

for sample in samples:
    sample_dictionary(sample)
    print(sample)

print("All done !")
