#!/usr/bin/env python3
# coding=utf-8

from itertools import chain
from collections import defaultdict
import numpy as np
from matplotlib import pyplot as plt
import os

### Interaction's dictionary
def dict_inter(data):
    interaction = {}
    with open("/home/laverre/Documents/Regulatory_Landscape/data/human/all_interactions/"+data, 'r') as f1:
        for i in f1.readlines()[1:]:
            i = i.split("\t")
            inter = (i[0] + "\t" + str(i[1]) + "\t" + str(i[2]) + "\t" + str(i[4]) + "\t" + str(i[5]) + "\t" + str(i[6]) + "\t")
            score = float(i[9].strip("\n"))
            #if i[0] == i[4]: only cis interactions
            interaction[inter] = score

    return interaction

dict1 = dict_inter("Bcell.ibed")
dict2 = dict_inter("CD34.ibed")
dict3 = dict_inter("PEK_early.ibed")
dict4 = dict_inter("PEK_late.ibed")
dict5 = dict_inter("PEK_undiff.ibed")
dict6 = dict_inter("pre_adipo.ibed")
dict7 = dict_inter("cardio.ibed")
dict8 = dict_inter("hESC.ibed")
dict9 = dict_inter("hNEC.ibed")
dict10 = dict_inter("MK.ibed")
dict11 = dict_inter("EP.ibed")
dict12 = dict_inter("Mon.ibed")
dict13 = dict_inter("TCD8.ibed")
dict14 = dict_inter("Ery.ibed")
dict15 = dict_inter("Neu.ibed")
dict16 = dict_inter("FoeT.ibed")
dict17 = dict_inter("NB.ibed")
dict18 = dict_inter("TCD4MF.ibed")
dict19 = dict_inter("TCD4Non.ibed")
dict20 = dict_inter("TCD4Act.ibed")
dict21 = dict_inter("Mac0.ibed")
dict22 = dict_inter("Mac1.ibed")
dict23 = dict_inter("Mac2.ibed")
dict24 = dict_inter("NCD4.ibed")
dict25 = dict_inter("TB.ibed")
dict26 = dict_inter("NCD8.ibed")



### Uniforming keys
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


### Fusion of interaction's dictionary
dict_final = defaultdict(list)
for k, v in chain(dict1.items(), dict2.items(), dict3.items(), dict4.items(), dict5.items(), dict6.items(),
                  dict7.items(), dict8.items(), dict9.items(), dict10.items(), dict11.items(), dict12.items(),
                  dict13.items(), dict14.items(), dict15.items(), dict16.items(), dict17.items(), dict18.items(),
                  dict19.items(), dict20.items(), dict21.items(), dict22.items(), dict23.items(), dict24.items(),
                  dict25.items(), dict26.items()):
    dict_final[k].append(v)


### Output
output = open("/home/laverre/Documents/Regulatory_Landscape/data/human/all_interactions/all_interactions.txt", 'w')
if os.stat("/home/laverre/Documents/Regulatory_Landscape/data/human/all_interactions/all_interactions.txt").st_size == 0:
    output.write("chr_bait\tstart_bait\tend_bait\tstart\tend\tBcell\tCD34\tPEK_early\tPEK_late\tPEK_undiff\tpre_adipo"
                 "\tcardio\thESC\thNEC\tMK\tEP\tMon\tTCD8\tEry\tNeu\tFoeT\tNB\tTCD4MF\tTCD4Non\tTCD4Act\tMac0\tMac1\t"
                 "Mac2\tNCD4\tTB\tNCD8\n")

for k, v in dict_final.items():
    output.write(k)
    for i in v:
        output.write(str(i) + "\t")
    output.write("\n")

output.close()


print("Nombre total d'interactions:", len(dict_final.keys()))
print("Nombre d'interaction unique:", (sum(len(i) == 1 for i in dict_final.values())))
print("Fréquence moyenne d'interaction :", (sum(len(i) for i in dict_final.values()))/len(dict_final.keys()))

freq = [len(i) for i in dict_final.values()]
print("Médiane fréquence d'interaction :", np.median(freq))

plt.hist(freq, min(len(set(freq)), 500), (0, 26))
plt.title('Human frequency')
plt.xlabel('cells')
plt.ylabel('Frequency')
plt.show()
plt.savefig('/home/laverre/Documents/Regulatory_Landscape/data/human/all_interactions/human_frequency.png')
