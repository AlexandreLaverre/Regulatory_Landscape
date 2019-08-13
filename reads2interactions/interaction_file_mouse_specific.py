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
    with open("/home/laverre/Documents/Regulatory_Landscape/data/mouse/all_interactions/"+data, 'r') as f1:
        for i in f1.readlines()[1:]:
            i = i.split("\t")
            inter = ('chr'+i[0] + "\t" + str(i[1]) + "\t" + str(i[2]) + "\t" + 'chr'+str(i[4]) + "\t" + str(i[5]) + "\t" + str(i[6]) + "\t")
            score = float(i[9].strip("\n"))
            #if i[0] == i[4]: only cis interactions
            interaction[inter] = score

    return interaction

dict2 = dict_inter("ESC.ibed")
dict3 = dict_inter("ESC_18.ibed")
dict4 = dict_inter("ESC_NKO.ibed")
dict5 = dict_inter("ESC_wild.ibed")


### Uniforming keys
L = [dict2, dict3, dict4, dict5]
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
for k, v in chain(dict2.items(), dict3.items(), dict4.items(), dict5.items()):
    dict_final[k].append(v)


### Output
output = open("/home/laverre/Documents/Regulatory_Landscape/data/mouse/all_interactions/all_interactions_chr_ESC.txt", 'w')
if os.stat("/home/laverre/Documents/Regulatory_Landscape/data/mouse/all_interactions/all_interactions_chr_ESC.txt").st_size == 0:
    output.write("chr_bait\tstart_bait\tend_bait\tstart\tend\tESC\tESC_18\tESC_NKO\tESC_wild\n")

for k, v in dict_final.items():
    output.write(k)
    for i in v:
        output.write(str(i) + "\t")
    output.write("\n")

output.close()

"""
print("Nombre total d'interactions:", len(dict_final.keys()))
print("Nombre d'interaction unique:", (sum(len(i) == 1 for i in dict_final.values())))
print("Fréquence moyenne d'interaction :", (sum(len(i) for i in dict_final.values()))/len(dict_final.keys()))

freq = [len(i) for i in dict_final.values()]
print("Médiane fréquence d'interaction :", np.median(freq))


### Frequency distribution of number of cells per interactions
plt.hist(freq, min(len(set(freq)), 500), (0, 15))
plt.title('Mouse frequency')
plt.xlabel('cells')
plt.ylabel('Frequency')
plt.show()
plt.savefig('/home/laverre/Documents/Regulatory_Landscape/data/mouse/all_interactions/mouse_frequency2.png')
"""
