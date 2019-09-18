#!/usr/bin/env python3
# coding=utf-8

import os
import numpy as np
# Conservation mouse interaction in human:
origin_sp = "human"  # or "human"
target_sp = "mouse"  # or "mouse"
data = "_simul"      # or "_simul"

frag_enh = {}
with open("../../data/"+origin_sp+"/overlap/"+origin_sp+"_frag_overlap_CAGE.txt") as f2:
    for i in f2.readlines()[1:]:
        i = i.strip("\n")
        i = i.split("\t")
        frag = ('chr' + str(i[0]) + ':' + str(i[1]) + ':' + str(i[2]))

        if i[3] != "NA":
            frag_enh[frag] = len(i[3].split(","))

cons_bait = {}
cons_bait_enh = {}
score_lift = {}
score_lift_enh = {}
print("Calculating align and synteny conservation... ")
with open("../../result/conservation/"+origin_sp+"2"+target_sp+"_conservation_syntenie_with_notconserv"+data+".txt2") as f1:
    for i in f1.readlines()[1:]:
        i = i.strip("\n")
        i = i.split("\t")
        bait = i[0].split('-')[0]
        PIR = i[0].split('-')[1]

        if bait not in cons_bait.keys():
            cons_bait[bait] = [0, 0, 0, 0, 0, 0, 0]  # (nb_PIR, nb_align, nb_synt, nb_inter, nb_TSS, nb_dvpt, nb_immun)
            cons_bait_enh[bait] = [0, 0, 0, 0]  # (nb_PIR, nb_align, nb_synt, nb_inter) enhancers !

        cons_bait[bait][0] += 1
        if PIR in frag_enh.keys():
            cons_bait_enh[bait][0] += 1

        if i[8] != 'NA':  # PIR lift
            cons_bait[bait][1] += 1
            if bait not in score_lift.keys():
                score_lift[bait] = [float(i[9])]
            else:
                score_lift[bait].append(float(i[9]))

            if PIR in frag_enh.keys():
                cons_bait_enh[bait][1] += 1
                if bait not in score_lift_enh.keys():
                    score_lift_enh[bait] = [float(i[9])]
                else:
                    score_lift_enh[bait].append(float(i[9]))

        if i[12] != 'NA':  # same chr
            cons_bait[bait][2] += 1
            if PIR in frag_enh.keys():
                cons_bait_enh[bait][2] += 1

med_lift = {}
med_lift_enh = {}
for bait in cons_bait.keys():
    if bait in score_lift.keys():
        med_lift[bait] = np.median(score_lift[bait])
    else:
        med_lift[bait] = 0

for bait in cons_bait.keys():
    if bait in score_lift_enh.keys():
        med_lift_enh[bait] = np.median(score_lift_enh[bait])
    else:
        med_lift_enh[bait] = 0

score_inter = {}
score_inter_enh = {}
print("Calculating interaction conservation... ")
with open("../../result/conservation/"+origin_sp+"2"+target_sp+"_conservation_interaction_pecan_0.4_merged"+data+".txt") as f2:
    for i in f2.readlines()[1:]:
        i = i.strip("\n")
        i = i.split("\t")
        bait = i[0].split('-')[0]
        PIR = i[0].split('-')[1]

        if i[5] != 'NA':  # interaction cons
            cons_bait[bait][3] += 1
            if PIR in frag_enh.keys():
                cons_bait_enh[bait][3] += 1
                if bait not in score_inter_enh.keys():
                    score_inter_enh[bait] = [float(i[11])]
                else:
                    score_inter_enh[bait].append(float(i[11]))

            if bait not in score_inter.keys():
                score_inter[bait] = [float(i[11])]
            else:
                score_inter[bait].append(float(i[11]))
            
med_inter = {}
for bait in cons_bait.keys():
    if bait in score_inter.keys():
        med_inter[bait] = np.median(score_inter[bait])
    else:
        med_inter[bait] = 0

med_inter_enh = {}
for bait in cons_bait.keys():
    if bait in score_inter_enh.keys():
        med_inter_enh[bait] = np.median(score_inter_enh[bait])
    else:
        med_inter_enh[bait] = 0

print("Calculating nb TSS involved in developmental & immun process... ")
dvpt = []
with open("../../data/"+origin_sp+"/annotations/GO_annotations/Gene_dvpt_process_QuickGO_uniq.txt") as f1:
    for i in f1.readlines():
        i = i.strip("\n")
        i = i.split("\t")
        dvpt.append(i[0])

immun = []
with open("../../data/"+origin_sp+"/annotations/GO_annotations/Gene_immune_process_QuickGO_uniq.txt") as f1:
    for i in f1.readlines():
        i = i.strip("\n")
        i = i.split("\t")
        immun.append(i[0])


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
                for x in TSS:
                    if x in dvpt:
                        cons_bait[frag][5] += 1
                    if x in immun:
                        cons_bait[frag][6] += 1

output = open("../../result/conservation/"+origin_sp+"2"+target_sp+"_conservation_by_bait_merged"+data+".txt", 'w')
if os.stat("../../result/conservation/"+origin_sp+"2"+target_sp+"_conservation_by_bait_merged"+data+".txt").st_size == 0:
    output.write("bait\tPIR\tPIR_enh\tseq_cons\tseq_cons_enh\tmed_cons\tmed_cons_enh\tPIR_synt\tPIR_synt_enh\t"
                 "PIR_int\tPIR_int_enh\tmed_inter\tmed_inter_enh\tTSS\tTSS_dvpt\tTSS_immun\n")

print("Writting output...")
for bait in cons_bait.keys():
    output.write(str(bait) + '\t' + str(cons_bait[bait][0]) + '\t' + str(cons_bait_enh[bait][0])
                 + '\t' + str(cons_bait[bait][1]) + '\t' + str(cons_bait_enh[bait][1])
                 + '\t' + str(med_lift[bait]) + '\t' + str(med_lift_enh[bait])
                 + '\t' + str(cons_bait[bait][2]) + '\t' + str(cons_bait_enh[bait][2])
                 + '\t' + str(cons_bait[bait][3]) + '\t' + str(cons_bait_enh[bait][3])
                 + '\t' + str(med_inter[bait]) + '\t' + str(med_inter_enh[bait])
                 + '\t' + str(cons_bait[bait][4]) + '\t' + str(cons_bait[bait][5]) + '\t' + str(cons_bait[bait][6]) + '\n')

output.close()
print("All done ! ")
