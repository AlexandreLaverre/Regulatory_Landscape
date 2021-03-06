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
            inter = ('chr' + i[0] + "\t" + str(i[1]) + "\t" + str(i[2]) + "\t" + 'chr'+str(i[4]) + "\t" + str(i[5]) + "\t" + str(i[6]) + "\t")
            score = float(i[9].strip("\n"))
            #if i[0] == i[4]: only cis interactions
            interaction[inter] = score

    return interaction


dict1 = dict_inter("Bcell.ibed")
dict2 = dict_inter("NB.ibed")
dict3 = dict_inter("TB.ibed")




### Uniforming keys
all = [dict1, dict2, dict3]

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
for k, v in chain(dict1.items(), dict2.items(), dict3.items()):
    dict_final[k].append(v)


### Output
output = open("/home/laverre/Documents/Regulatory_Landscape/data/human/all_interactions/all_interactions_chr_Bcell.txt", 'w')
if os.stat("/home/laverre/Documents/Regulatory_Landscape/data/human/all_interactions/all_interactions_chr_Bcell.txt").st_size == 0:
    output.write("chr_bait\tstart_bait\tend_bait\tstart\tend\tBcell\tNB\tTB\n")

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

plt.hist(freq, min(len(set(freq)), 500), (0, 26))
plt.title('Human frequency')
plt.xlabel('cells')
plt.ylabel('Frequency')
plt.show()
plt.savefig('/home/laverre/Documents/Regulatory_Landscape/data/human/all_interactions/human_frequency2.png')
"""