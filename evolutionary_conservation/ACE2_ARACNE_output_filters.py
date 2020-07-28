#!/usr/bin/env python3
# coding=utf-8

import os

input_file = "/home/laverre/Documents/ACE2/results/mutual_information_network/navratil/ACE2_ARACNE_results.txt"
output_file = "/home/laverre/Documents/ACE2/results/mutual_information_network/navratil/ACE2_ARACNE_results_summary2.csv"

ID_to_name = {}
with open("/home/laverre/Documents/ACE2/data/ensembl_annotations/GeneNames_Ensembl99.txt") as f1:
    for i in f1.readlines()[1:]:
        i = i.strip("\n")
        i = i.split("\t")
        ID_to_name[i[0]] = str(i[1])

ACE2_TF = {"test": {"test": "test"}}
with open(input_file) as f3:
    for i in f3.readlines():
        i = i.strip("\n")
        i = i.split("\t")
        group = i[0].split("/")[0]
        TF = ID_to_name[i[2]]
        MI = i[4]
        pearson = i[5]

        if "nan" not in i: #pval
            if TF in ACE2_TF.keys():
                ACE2_TF[TF][str(group)] = str(MI + ';' + pearson)
            else:
                ACE2_TF[TF] = {}
                ACE2_TF[TF][str(group)] = str(MI + ';' + pearson)

groups = ['group_A', 'group_B', 'group_C', 'group_D', 'group_E', 'group_F', 'group_G', 'group_H', 'group_I', 'group_J', 'group_K']

output = open(output_file, 'w')
if os.stat(output_file).st_size == 0:
    output.write("TF" + "\t" + '\t'.join(groups) + "\n")


for TF, TF_group in ACE2_TF.items():
    output.write(TF + "\t")
    for group in groups:
        if group == "group_K":
            if group in TF_group:
                output.write(TF_group[group] + '\n')
            else:
                output.write("" + '\n')

        else:
            if group in TF_group:
                output.write(TF_group[group] + '\t')
            else:
                output.write("" + '\t')

output.close()

