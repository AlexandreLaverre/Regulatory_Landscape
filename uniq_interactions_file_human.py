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
    with open("/home/laverre/Documents/Regulatory_Landscape/data/human/all_interactions/ibed/"+data, 'r') as f1:
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


### Uniforming keys
L = [dict1, dict2, dict3, dict4, dict5, dict6, dict7, dict8, dict9]
alldicts = {}
for d in L:
    alldicts.update(d)
allkeys = alldicts.keys()

for d in L:
    for key in allkeys:
        if key not in d:
            d[key] = 'NA'

### Fusion of interaction's dictionary
dict_final = defaultdict(list)
for k, v in chain(dict1.items(), dict2.items(), dict3.items(), dict4.items(), dict5.items(), dict6.items(), dict7.items(), dict8.items(), dict9.items()):
    dict_final[k].append(v)

### Output
output = open("/home/laverre/Documents/Regulatory_Landscape/data/human/all_interactions/ibed/all_interactions.txt", 'w')
if os.stat("/home/laverre/Documents/Regulatory_Landscape/data/human/all_interactions/ibed/all_interactions.txt").st_size == 0:
    output.write("chr_bait\tstart_bait\tend_bait\tstart\tend\tBcell\tCD34\tPEK_early\tPEK_late\tPEK_undiff\tpre_adipo\tcardio\thESC\thNEC\n")

for k, v in dict_final.items():
    output.write(k)
    for i in v:
        output.write(str(i) + "\t")
    output.write("\n")

print("Nombre total d'interactions:", len(dict_final.keys()))
print("Nombre d'interaction unique:", (sum(len(i) == 1 for i in dict_final.values())))
print("Fréquence moyenne d'interaction :", (sum(len(i) for i in dict_final.values()))/len(dict_final.keys()))

freq = [len(i) for i in dict_final.values()]
print("Médiane fréquence d'interaction :", np.median(freq))

plt.hist(freq, min(len(set(freq)), 500), (0, 10))
plt.title('Human frequency')
plt.xlabel('cells')
plt.ylabel('Frequency')
plt.show()
