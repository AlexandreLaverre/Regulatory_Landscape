#!/usr/bin/env python3
# coding=utf-8

import os

# Conservation mouse interaction in human:
origin_sp = "human"
target_sp = "mouse"

cons_bait = {}
print("Calculating align and synteny conservation... ")
with open("../../result/conservation/"+origin_sp+"2"+target_sp+"_conservation_syntenie_with_notconserv.txt") as f1:
    for i in f1.readlines()[1:]:
        i = i.strip("\n")
        i = i.split("\t")
        bait = i[0].split('-')[0]

        if bait not in cons_bait.keys():
            cons_bait[bait] = [0, 0, 0, 0, 0]  # (nb_PIR, nb_align, nb_synt, nb_inter)

        cons_bait[bait][0] += 1
        if i[6] != 'NA':  # PIR lift
            cons_bait[bait][1] += 1

        if i[8] != 'NA':  # same chr
            cons_bait[bait][2] += 1

print("Calculating interaction conservation... ")
with open("../../result/conservation/"+origin_sp+"2"+target_sp+"_conservation_interaction.txt") as f2:
    for i in f2.readlines()[1:]:
        i = i.strip("\n")
        i = i.split("\t")
        bait = i[0].split('-')[0]

        if i[5] != 'NA':  # interaction cons
            cons_bait[bait][3] += 1

print("Calculating nb TSS... ")
with open("../../data/"+origin_sp+"/overlap/"+origin_sp+"_frag_overlap_TSS.txt") as f1:
    for i in f1.readlines()[1:]:
        i = i.strip("\n")
        i = i.split("\t")
        frag = "chr"+i[0]+':'+i[1]+':'+i[2]

        if frag in cons_bait.keys():
            TSS = i[3].split(',')

            if TSS != ['NA']:
                cons_bait[frag][4] = len(TSS)


output = open("../../result/conservation/"+origin_sp+"2"+target_sp+"_conservation_by_bait.txt", 'w')
if os.stat("../../result/conservation/"+origin_sp+"2"+target_sp+"_conservation_by_bait.txt").st_size == 0:
    output.write("bait\tPIR\tPIR_lift\tPIR_synt\tPIR_int\tTSS\n")

print("Writting output...")
for bait in cons_bait.keys():
    output.write(str(bait) + '\t' + str(cons_bait[bait][0]) + '\t' + str(cons_bait[bait][1]) + '\t'
                 + str(cons_bait[bait][2]) + '\t' + str(cons_bait[bait][3]) + '\t' + str(cons_bait[bait][4]) + '\n')

output.close()
print("All done ! ")
