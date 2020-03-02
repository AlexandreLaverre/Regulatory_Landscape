#!/usr/bin/env python3
# coding=utf-8

from itertools import chain
from collections import defaultdict
import numpy as np
from matplotlib import pyplot as plt
import os


origin = "simulated_"
sp = "mouse_"
path = "/home/laverre/Documents/Regulatory_Landscape/data/"

if origin == "observed":
    data = ""
    extension = ".ibed"
    path_data = path + "/mouse/all_interactions/mouse_samples/"

else:
    data = "simulated_"
    extension = "_bait_other.txt"
    path_data = path + "/mouse/Simulations/simulations_samples/"

# Baited fragment
baits = {}
with open(path + "/mouse/Digest_HindIII_None.txt.baitmap") as f1:
    for i in f1.readlines():
        i = i.strip("\n")
        i = i.split("\t")
        frag = "chr" + str(i[0]) + "\t" + str(i[1]) + "\t" + str(i[2])
        baits[frag] = "bait"


### Interaction's dictionary
def dict_inter(sample):
    interaction = {}
    with open(path_data+data+sample+extension, 'r') as f1:
        for i in f1.readlines()[1:]:
            i = i.strip("\n")
            i = i.split("\t")

            bait = "chr" + str(i[0]) + "\t" + str(i[1]) + "\t" + str(i[2])
            midbait = (int(i[1]) + int(i[2])) / 2
            contacted = "chr" + str(i[3]) + "\t" + str(i[4]) + "\t" + str(i[5])
            midcontacted = (int(i[4]) + int(i[5])) / 2

            dist = str(abs(midbait - midcontacted))
            type = "baited" if contacted in baits.keys() else "unbaited"

            contact = (bait + "\t" + contacted + "\t" + type + "\t" + dist + "\t")

            if origin == "observed":
                score = float(i[7].strip("\n"))
            else:
                score = 1

            interaction[contact] = score

    return interaction


dict1 = dict_inter("EpiSC")
dict2 = dict_inter("ESC")
dict3 = dict_inter("ESC_18")
dict4 = dict_inter("ESC_NKO")
dict5 = dict_inter("ESC_wild")
dict6 = dict_inter("ESd_starved")
dict7 = dict_inter("ESd_TPO")
dict8 = dict_inter("FLC")
dict9 = dict_inter("preB_aged")
dict10 = dict_inter("preB_young")
dict11 = dict_inter("TSC")
dict12 = dict_inter("preadip_D0")
dict13 = dict_inter("preadip_D2")
dict14 = dict_inter("preadip_4H")

### Uniforming keys
all = [dict1, dict2, dict3, dict4, dict5, dict6, dict7, dict8, dict9, dict10, dict11, dict12, dict13, dict14]
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
                  dict13.items(), dict14.items()):
    dict_final[k].append(v)


### Output
output = open(path_data+origin+"_all_interactions.txt", 'w')
if os.stat(path_data+origin+"_all_interactions.txt").st_size == 0:
    output.write("chr_bait\tstart_bait\tend_bait\tchr\tstart\tend\ttype\tdistance\tEpiSC\tESC\tESC_18\tESC_NKO\tESC_wild\tESd_starved\tESd_TPO"
                 "\tFLC\tpreB_aged\tpreB_young\tTSC\tpreadip_D0\tpreadip_D2\tpreadip_4H\n")

for k, v in dict_final.items():
    output.write(k)
    for i in v:
        output.write(str(i) + "\t")
    output.write("\n")

output.close()

print("Nombre total d'interactions:", len(dict_final.keys()))
print("Nombre d'interaction unique:", (sum(i.count("NA") == len(i)-1 for i in dict_final.values())))
print("Fréquence moyenne d'interaction :", (sum(i.count(1) for i in dict_final.values()))/len(dict_final.keys()))

freq = [len(i) for i in dict_final.values()]
print("Médiane fréquence d'interaction :", np.median(freq))


# ### Frequency distribution of number of cells per interactions
# plt.hist(freq, min(len(set(freq)), 500), (0, 15))
# plt.title('Mouse frequency')
# plt.xlabel('cells')
# plt.ylabel('Frequency')
# plt.show()
# plt.savefig('/home/laverre/Documents/Regulatory_Landscape/data/mouse/all_interactions/mouse_frequency2.png')

