#!/usr/bin/env python3
# coding=utf-8

import os

# Conservation mouse interaction in human:
origin_sp = "mouse"  # or "human"
target_sp = "human"  # or "mouse"
data = ""  # or "_simul"


print("Calculating align and synteny conservation... ")

cons_bait = {}
with open("../../result/conservation/" + origin_sp + "2" + target_sp + "_conservation_syntenie_with_notconserv"
          + data + ".txt") as f1:
    for i in f1.readlines()[1:]:
        i = i.strip("\n")
        i = i.split("\t")
        bait = i[0].split('-')[0]

        if bait not in cons_bait.keys():
            cons_bait[bait] = [0, 0, 0, 0]  # (nb_PIR, nb_align, nb_synt, nb_inter)

        cons_bait[bait][0] += 1
        if i[8] != 'NA':  # PIR lift
            cons_bait[bait][1] += 1

        if i[12] != 'NA':  # PIR in synteny
            cons_bait[bait][2] += 1


print("Calculating interaction conservation... ")

with open("../../result/conservation/" + origin_sp + "2" + target_sp + "_conservation_interaction" + data + ".txt") as f2:
    for i in f2.readlines()[1:]:
        i = i.strip("\n")
        i = i.split("\t")
        bait = i[0].split('-')[0]

        if i[5] != 'NA':  # interaction cons
            cons_bait[bait][3] += 1


print("Linking bait to TSS gene ID... ")

cons_gene = {}
with open("../../data/" + origin_sp + "/overlap/" + origin_sp + "_frag_overlap_TSS.txt2") as f1:
    for i in f1.readlines()[1:]:
        i = i.strip("\n")
        i = i.split("\t")
        frag = "chr" + i[0] + ':' + i[1] + ':' + i[2]

        if frag in cons_bait.keys():  # select only baited fragment
            all_TSS = i[3].split(',')

            if all_TSS != ['NA']:       # select only baited fragment overlap with TSS
                for TSS in all_TSS:
                    if TSS not in cons_gene.keys():
                        cons_gene[TSS] = [cons_bait[frag]]
                    else:
                        cons_gene[TSS].append(cons_bait[frag])

sum_cons_gene = {}
for TSS in cons_gene.keys():
    if len(cons_gene[TSS]) > 1:
        sum_cons_gene[TSS] = [sum(x) for x in zip(*cons_gene[TSS])]

    else:
        sum_cons_gene[TSS] = cons_gene[TSS][0]


print("Writting output...")

output = open("../../result/conservation/" + origin_sp + "2" + target_sp + "_conservation_by_gene" + data + ".txt", 'w')
if os.stat("../../result/conservation/" + origin_sp + "2" + target_sp + "_conservation_by_gene" + data + ".txt").st_size == 0:
    output.write("TSS\tnb_PIR\tnb_PIR_lift\tnb_PIR_synt\tnb_PIR_int\n")

for TSS in sum_cons_gene.keys():
    output.write(str(TSS) + '\t' + str(sum_cons_gene[TSS][0]) + '\t' + str(sum_cons_gene[TSS][1]) + '\t' +
                 str(sum_cons_gene[TSS][2]) + '\t' + str(sum_cons_gene[TSS][3]) + '\n')

output.close()
print("All done ! ")
