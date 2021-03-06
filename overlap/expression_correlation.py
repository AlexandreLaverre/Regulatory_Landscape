#!/usr/bin/env python3
# coding=utf-8

from scipy import stats
import numpy as np
import os
from datetime import datetime
import sys

startTime = datetime.now()

path_data = "/home/laverre/Documents/Regulatory_Landscape/data/human/"
path_overlap = "/home/laverre/Documents/Regulatory_Landscape/data/human/overlap/"
path_result = "/home/laverre/Documents/Regulatory_Landscape/result/expression_correlation/"

"""
path_data = "/beegfs/data/alaverre/data/expression_correlation/"+sp+"/"
path_overlap = "/beegfs/data/alaverre/data/fragments_overlap/"+sp+"/"
path_result = "/beegfs/data/alaverre/result/expression_correlation/"+sp+"/"
"""

def frag_dictionary(file):
    dic = {}
    with open(file, 'r') as f:
        for i in f.readlines()[1:]:
            i = i.strip("\n")
            i = i.split("\t")
            frag = (i[0] + ":" + i[1] + ":" + i[2])  # key = fragment coord
            ID = i[3].split(",")                     # value = overlap_ID

            if len(ID) == 1:
                dic[frag] = ID[0]

            else:
                for x in ID:
                    if frag in dic.keys():
                        dic[frag].append(x)
                    else:
                        dic[frag] = [x]
    return dic


def expr_dictionary(file, first_cell):
    dic = {}
    with open(file, 'r') as f:
        for i in f.readlines()[1:]:
            i = i.strip('\n')
            i = i.split('\t')
            ID = str(i[0])                          # key = peaks ID

            dic[ID] = []
            for x in i[first_cell:len(i)]:
                dic[ID].append(float(x))            # value = peaks expression
    return dic


# Restriction fragment - CAGE peaks
frag_peaks = frag_dictionary(path_overlap+"human_overlap_CAGE_peaks_filtered.txt")
print("Restriction fragment - CAGE peaks ok !")

# Restriction fragment - enhancers peaks
frag_enh = frag_dictionary(path_overlap+"human_overlap_CAGE_enh_filtered.txt")
print("Restriction fragment - enhancers peaks ok !")

# CAGE peaks expression
cell = 1
exp_CAGE = expr_dictionary(path_data+"/CAGE/human.CAGE_peaks.expression.tpm.ann.matrix_filtered", cell)
print("CAGE peaks expression ok !")

# enhancers expression
cell = 1
exp_enh = expr_dictionary(path_data+"/CAGE/human.enhancers.expression.tpm.matrix_filtered", cell)
print("Enhancers peaks expression ok !")

# pc-HIC interactions
inter = {}
with open(path_data+"/Simulations/simulations_human_10Mb_bin5kb_fragoverbin_chr.txt", 'r') as f:
    for i in f.readlines()[1:]:
        i = i.strip("\n")
        i = i.split("\t")
        bait = (i[0] + ":" + i[1] + ":" + i[2])  # keys = bait
        midbait = ((int(i[1]) + int(i[2])) / 2)
        PIR = (i[3] + ":" + i[4] + ":" + i[5])   # value = PIR
        midPIR = ((int(i[4]) + int(i[5])) / 2)

        if i[0] == i[3]:
            if 25000 < abs(midbait-midPIR) < 10000000:
                if bait in inter.keys():
                    inter[bait].append(PIR)
                else:
                    inter[bait] = [PIR]

# Correlation

print("Making correlations... ")

output = open(path_result+"expression_correlation.txt_new_simul", 'w')
if os.stat(path_result+"expression_correlation.txt_new_simul").st_size == 0:
    output.write("peaks\tenh\tspearman_cor\tspearman_pval\tpearson_cor\tpearson_pval\tbait\tPIR\tsense\n")


nb_pairs = 0
for bait in inter.keys():
    if frag_peaks[bait] != "NA":        # Get CAGE peaks in bait
        for PIR in inter[bait]:
            if frag_enh[PIR] != "NA":       # Get enh peaks in PIR
                for peak in frag_peaks[bait]:
                    if peak in exp_CAGE.keys():
                        vect_peak = exp_CAGE[peak]  # Get expression of CAGE peak
                        for enh in frag_enh[PIR]:
                            if enh in exp_enh.keys():
                                vect_enh = exp_enh[enh]  # Get expression of enh peak

                                cor_spearman = stats.spearmanr(vect_peak, vect_enh)
                                cor_pearson = stats.pearsonr([np.log(i+1) for i in vect_peak],
                                                             [np.log(i+1) for i in vect_enh])

                                nb_pairs += 1
                                output.write(peak + "\t" + enh + "\t" + str(cor_spearman[0]) + "\t" +
                                             str(cor_spearman[1]) + "\t" + str(cor_pearson[0]) + "\t" +
                                             str(cor_pearson[1]) + "\t" + bait + "\t" + PIR + "\t" + "bait-PIR" + "\n")

    if frag_enh[bait] != "NA":      # Check enh peaks in bait
        for PIR in inter[bait]:
            if frag_peaks[PIR] != "NA":     # Check CAGE peaks in PIR
                for enh in frag_enh[bait]:
                    if enh in exp_enh.keys():
                        vect_enh = exp_enh[enh]  # Get expression of enh peak
                        for peak in frag_peaks[PIR]:
                            if peak in exp_CAGE.keys():
                                vect_peak = exp_CAGE[peak]   # Get expression of CAGE peak
                                cor_spearman = stats.spearmanr(vect_peak, vect_enh)
                                cor_pearson = stats.pearsonr([np.log(i + 1) for i in vect_peak],
                                                             [np.log(i + 1) for i in vect_enh])

                                nb_pairs += 1
                                output.write(peak + "\t" + enh + "\t" + str(cor_spearman[0]) + "\t" +
                                             str(cor_spearman[1]) + "\t" + str(cor_pearson[0]) + "\t" +
                                             str(cor_pearson[1]) + "\t" + bait + "\t" + PIR + "\t" + "PIR-bait" + "\n")


print("Pairs number:", nb_pairs)
print("All done !")
print("Execution time :", datetime.now() - startTime)
output.close()

