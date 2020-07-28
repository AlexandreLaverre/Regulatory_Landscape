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
    all_interactions = path_data + "/Simulations/simulations_" + sp + "_10Mb_bin5kb_fragoverbin_chr_info.txt"
    #all_interactions = path_data + "simulated_all_interactions.txt"
else:
    data_type = "observed"
    all_interactions = path_data + "/all_interactions/all_interactions_chr_info.txt"

if sp == "mouse":
    genome = "mm10"
else:
    genome = "hg38"

print("specie:", sp)
print("data:", data_type)

# PIR overlap bait
PIR_bait = {}
with open(
        "/home/laverre/Data/Regulatory_landscape/data/CHICAGO_files/" + sp + "/" + genome + "/Digest_" + genome + "_HindIII_None.txt.baitmap") as f1:
    for i in f1.readlines():
        i = i.strip("\n")
        i = i.split("\t")
        frag = "chr" + str(i[0]) + "\t" + str(i[1]) + "\t" + str(i[2])
        PIR_bait[frag] = "bait"


# Sample interaction dic
def sample_dictionary(sample):
    sample_dic = {}
    with open(sample, 'r') as f:
        for i in f.readlines()[1:]:
            i = i.strip("\n")
            i = i.split("\t")
            bait = ("chr" + str(i[0]) + "\t" + str(i[1]) + "\t" + str(i[2]))

            if data == "_simul":
                PIR = ("chr" + str(i[3]) + "\t" + str(i[4]) + "\t" + str(i[5]), "NA", "NA")
            else:
                PIR = ("chr" + str(i[3]) + "\t" + str(i[4]) + "\t" + str(i[5]), str(i[6]),
                       str(i[7]))  # chr, start, end, N_reads, score

            if bait not in sample_dic.keys():
                sample_dic[bait] = [PIR]
            else:
                sample_dic[bait].append(PIR)

    # Writting output
    output_file = sample + "_infos"
    output = open(output_file, 'w')
    if os.stat(output_file).st_size == 0:
        output.write(
            "bait_chr\tbait_start\tbait_end\tchr\tstart\tend\tN_reads\tscore\tbaited_frag\tdist\n")

    for bait, contacted in sample_dic.items():
        for i in contacted:
            if bait.split('\t')[0] == i[0].split('\t')[0]:
                midbait = ((int(bait.split('\t')[1]) + int(bait.split('\t')[2])) / 2)
                midcontact = ((int(i[0].split('\t')[1]) + int(i[0].split('\t')[2])) / 2)
                print(midbait, midcontact)
                dist_obs = abs(midbait - midcontact)
                print(dist_obs)
            else:
                dist_obs = "NA"

            PIR = str(i[0])
            PIR_status = "baited" if PIR in PIR_bait.keys() else "unbaited"

            output.write(
                bait + "\t" + PIR + "\t" + i[1] + "\t" + i[2] + "\t" + PIR_status + "\t" + str(dist_obs) + "\n")

    output.close()


if data == "_simul":
    samples = glob.glob(path_data + "Simulations/simulations_samples/*txt")
else:
    samples = glob.glob(path_data + "all_interactions/" + sp + "_samples/*ibed")

for sample in samples:
    sample_dictionary(sample)
    print(sample)

print("All done !")
