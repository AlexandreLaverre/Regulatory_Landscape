#!/usr/bin/env python3
# coding=utf-8

import os

# Conservation mouse interaction in human:
origin_sp = "mouse"  # or "human"
target_sp = "human"  # or "mouse"
data = ""  # or "_simul"
cell_type = "Bcell"


dic_specific = {}
with open("../../data/"+origin_sp+"/all_interactions/all_interactions_chr_"+cell_type+".txt") as specific:
    for i in specific.readlines()[1:]:
        i = i.strip("\n")
        i = i.split("\t")
        interaction = str(i[0]+':'+i[1]+':'+i[2]+'-'+i[3]+':'+i[4]+':'+i[5])
        dic_specific[interaction] = 0


print("Calculating align and synteny conservation... ")

cons_bait = {}
with open("../../result/conservation/" + origin_sp + "2" + target_sp + "_conservation_syntenie_with_notconserv"
          + data + ".txt2") as f1:
    for i in f1.readlines()[1:]:
        i = i.strip("\n")
        i = i.split("\t")
        interaction = i[0]

        if interaction in dic_specific.keys():
            bait = i[0].split('-')[0]

            if bait not in cons_bait.keys():
                cons_bait[bait] = [0, 0, 0, 0]  # (nb_PIR, nb_align, nb_synt, nb_inter)

            cons_bait[bait][0] += 1  # total PIR

            if i[9] != 'NA':
                if float(i[9]) >= 0.4:  # conserv seq PIR
                    cons_bait[bait][1] += 1

            if i[12] != 'NA':
                if float(i[12]) < 10000000:  # conserv synteny PIR
                    cons_bait[bait][2] += 1


print("Calculating interaction conservation... ")

with open("../../result/conservation/" + origin_sp + "2" + target_sp + "_conservation_interaction_pecan_0.4_" + cell_type + data + ".txt") as f2:
    for i in f2.readlines()[1:]:
        i = i.strip("\n")
        i = i.split("\t")
        interaction = i[0]

        if interaction in dic_specific.keys():
            bait = i[0].split('-')[0]

            if i[5] != 'NA':  # interaction cons
                cons_bait[bait][3] += 1


print("Linking bait to TSS gene ID... ")

cons_gene = {}
with open("../../data/" + origin_sp + "/overlap/" + origin_sp + "_frag_overlap_TSS.txt") as f1:
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

output = open("../../result/conservation/" + origin_sp + "2" + target_sp + "_conservation_by_gene_0.4_" + cell_type + data + ".txt2", 'w')
if os.stat("../../result/conservation/" + origin_sp + "2" + target_sp + "_conservation_by_gene_0.4_" + cell_type + data + ".txt2").st_size == 0:
    output.write("TSS\tnb_contact\tnb_seq_conserv\tnb_synt_conserv\tnb_int_conserv\n")

for TSS in sum_cons_gene.keys():
    output.write(str(TSS) + '\t' + str(sum_cons_gene[TSS][0]) + '\t' + str(sum_cons_gene[TSS][1]) + '\t' +
                 str(sum_cons_gene[TSS][2]) + '\t' + str(sum_cons_gene[TSS][3]) + '\n')

output.close()
print("All done ! ")
