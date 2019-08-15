#!/usr/bin/env python3
# coding=utf-8

from itertools import chain
from collections import defaultdict
import os


### Interaction's dictionary
def dict_inter(data):
    dic = {}
    with open("/home/laverre/Documents/Regulatory_Landscape/data/human/all_interactions/"+data, 'r') as f1:
        for i in f1.readlines()[1:]:
            i = i.strip("\n")
            i = i.split("\t")

            if i[0] == i[4]:  # only cis interactions
                bait = (str(i[0]) + "\t" + str(i[1]) + "\t" + str(i[2]))
                score = float(i[9].strip("\n"))
                PIR = (i[4], i[5], i[6], score)

                midbait = ((int(i[1]) + int(i[2])) / 2)
                midcontact = ((int(i[5]) + int(i[6])) / 2)
                dist = abs(midbait - midcontact)

                if 25000 < dist <= 10000000:
                    if bait in dic.keys():
                        dic[bait].append(PIR)
                    else:
                        dic[bait] = [PIR]

        # Sorting tuples for each key
        for k in dic.keys():
            dic[k] = list(set(tuple(x) for x in dic[k]))
            dic[k].sort(key=lambda x: x[1])

        # Merging adjacent contacted fragments
        merged_dic = {}
        for bait in dic.keys():
            chr = dic[bait][0][0]
            current_start = dic[bait][0][1]
            current_end = dic[bait][0][2]
            current_score = dic[bait][0][3]
            merged_dic[bait] = []
            merge_info = "no_merged"
            merge_count = 1

            if len(dic[bait]) == 1:
                merged_dic[bait].append((chr, current_start, current_end, current_score, merge_info))
            else:
                for i in range(1, len(dic[bait])):
                    new_start = dic[bait][i][1]
                    new_end = dic[bait][i][2]
                    new_score = dic[bait][i][3]

                    if int(new_start) == int(current_end) + 1:
                        merge_count += 1
                        current_end = new_end
                        merge_info = "merged"
                        current_score = current_score+new_score

                    else:
                        current_score = current_score / merge_count
                        merged_dic[bait].append((chr, current_start, current_end, current_score, merge_info))
                        current_start = new_start
                        current_end = new_end
                        merge_info = "no_merged"
                        merge_count = 1

                    if i == len(dic[bait]) - 1:
                        merged_dic[bait].append((chr, current_start, current_end, current_score, merge_info))

        interactions = {}
        for bait in merged_dic.keys():
            for contact in merged_dic[bait]:
                inter = (bait + "\t" + str(contact[0]) + "\t" + str(contact[1]) + "\t" + str(contact[2]) + "\t" + str(contact[4]) + "\t")
                score = str(contact[3])
                interactions[inter] = score

        return interactions


dict1 = dict_inter("Bcell.ibed")
dict2 = dict_inter("CD34.ibed")
dict3 = dict_inter("PEK_early.ibed")
dict4 = dict_inter("PEK_late.ibed")
dict5 = dict_inter("PEK_undiff.ibed")
print("Merging sample : 5 / 26")
dict6 = dict_inter("pre_adipo.ibed")
dict7 = dict_inter("cardio.ibed")
dict8 = dict_inter("hESC.ibed")
dict9 = dict_inter("hNEC.ibed")
dict10 = dict_inter("MK.ibed")
print("Merging sample : 10 / 26")
dict11 = dict_inter("EP.ibed")
dict12 = dict_inter("Mon.ibed")
dict13 = dict_inter("TCD8.ibed")
dict14 = dict_inter("Ery.ibed")
dict15 = dict_inter("Neu.ibed")
print("Merging sample : 15 / 26")
dict16 = dict_inter("FoeT.ibed")
dict17 = dict_inter("NB.ibed")
dict18 = dict_inter("TCD4MF.ibed")
dict19 = dict_inter("TCD4Non.ibed")
dict20 = dict_inter("TCD4Act.ibed")
print("Merging sample : 20 / 26")
dict21 = dict_inter("Mac0.ibed")
dict22 = dict_inter("Mac1.ibed")
dict23 = dict_inter("Mac2.ibed")
dict24 = dict_inter("NCD4.ibed")
dict25 = dict_inter("TB.ibed")
dict26 = dict_inter("NCD8.ibed")
print("Merging sample done ! ")

# Uniforming keys
all = [dict1, dict2, dict3, dict4, dict5, dict6, dict7, dict8, dict9, dict10, dict11, dict12, dict13, dict14, dict15,
     dict16, dict17, dict18, dict19, dict20, dict21, dict22, dict23, dict24, dict25, dict26]

alldicts = {}
for dic in all:
    alldicts.update(dic)
allkeys = alldicts.keys()

for dic in all:
    for key in allkeys:
        if key not in dic:
            dic[key] = 'NA'

# Fusion of interaction's dictionary
dict_final = defaultdict(list)
for k, v in chain(dict1.items(), dict2.items(), dict3.items(), dict4.items(), dict5.items(), dict6.items(),
                  dict7.items(), dict8.items(), dict9.items(), dict10.items(), dict11.items(), dict12.items(),
                  dict13.items(), dict14.items(), dict15.items(), dict16.items(), dict17.items(), dict18.items(),
                  dict19.items(), dict20.items(), dict21.items(), dict22.items(), dict23.items(), dict24.items(),
                  dict25.items(), dict26.items()):
    dict_final[k].append(v)

# Output
print("Writting output...")
output = open("/home/laverre/Documents/Regulatory_Landscape/data/human/all_interactions/all_interactions_merged.txt", 'w')
if os.stat("/home/laverre/Documents/Regulatory_Landscape/data/human/all_interactions/all_interactions_merged.txt").st_size == 0:
    output.write("chr_bait\tstart_bait\tend_bait\tchr\tstart\tend\tmerge_info\tBcell\tCD34\tPEK_early\tPEK_late\tPEK_undiff\tpre_adipo"
                 "\tcardio\thESC\thNEC\tMK\tEP\tMon\tTCD8\tEry\tNeu\tFoeT\tNB\tTCD4MF\tTCD4Non\tTCD4Act\tMac0\tMac1\t"
                 "Mac2\tNCD4\tTB\tNCD8\n")

for k, v in dict_final.items():
    output.write(k)
    for i in v:
        output.write(str(i) + "\t")
    output.write("\n")

output.close()
