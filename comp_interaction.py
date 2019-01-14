#!/usr/bin/env python3
# coding=utf-8

# Dictionnaire interactions Cairns
bait_PIR = {}
tot_inter = 0

with open("GSE81503_mESC_PCHiC_merge_final_washU_text2.txt", 'r') as f1:
    for i in f1.readlines()[1:]:
        i = i.strip('\n')
        i = i.split("\t")
        bait = i[0]
        PIR = i[1]
        tot_inter += 1
        if bait in bait_PIR.keys():
            bait_PIR[bait].append(PIR)
        if PIR in bait_PIR.keys():
            bait_PIR[PIR].append(bait)
        else:
            bait_PIR[bait] = [PIR]
            bait_PIR[PIR] = [bait]

print("Interactions totales in Cairns :", sum(len(bait_PIR[i]) for i in bait_PIR.keys()))


# Dictionnaire interactions perso
bait_PIR_perso = {}
PIR_bait_perso = {}
tot_inter_perso = 0

with open("ESC_sorted_dedup_washU_text.txt", 'r') as f1:
    for i in f1.readlines()[1:]:
        i = i.strip('\n')
        i = i.split("\t")
        bait = i[0]
        PIR = i[1]
        tot_inter_perso += 1
        if bait in bait_PIR_perso.keys():
            bait_PIR_perso[bait].append(PIR)
        if PIR in bait_PIR_perso.keys():
            bait_PIR_perso[PIR].append(bait)
        else:
            bait_PIR_perso[bait] = [PIR]
            bait_PIR_perso[PIR] = [bait]

print("Interactions totales perso :", sum(len(bait_PIR_perso[i]) for i in bait_PIR_perso.keys()))

# Comparaison
interaction_ok = 0
interaction_non_ok = 0

for bait in bait_PIR.keys():
    if bait in bait_PIR_perso.keys():
        for PIR in bait_PIR[bait]:
            if PIR in bait_PIR_perso[bait]:
                interaction_ok += 1
            else:
                interaction_non_ok += 1
    else:
        for PIR in bait_PIR[bait]:
            interaction_non_ok += 1



print("Interactions communes:", interaction_ok)
print("Autres:", interaction_non_ok)
