#!/usr/bin/env python3
# coding=utf-8

import os


def interaction_dictionary(file):
    inter = {}
    dist = []
    with open(file) as f3:
        for i in f3.readlines()[1:]:
            i = i.strip("\n")
            i = i.split("\t")
            frag = (str(i[0]) + ":" + str(i[1]) + "-" + str(i[2]))
            PIR = (str(i[3]) + ':' + str(i[4]) + "-" + str(i[5]))
            midbait = ((int(i[1]) + int(i[2])) / 2)
            midcontact = ((int(i[4]) + int(i[5])) / 2)

            dist.append(midbait - midcontact)

            if frag in inter.keys():
                inter[frag].append(PIR)
            else:
                inter[frag] = [PIR]

    return inter, dist


# Dict lifted interaction
inter_lift, dist_lift = interaction_dictionary("../../result/alignments/mouse2human/mouse2human_all_interaction_simul.txt")

# Dict real interaction
inter_real, dist_real = interaction_dictionary("../../data/human/all_interactions/all_interactions_chr.txt")

# Dict lifted_frag == real_frag
lift_real = {}
with open("../../data/human/overlap/mouse2human_frag_overlap_human_frag.txt") as f2:
    for i in f2.readlines()[1:]:
        i = i.strip("\n")
        i = i.split("\t")
        lifted_frag = (str(i[0]) + ":" + str(i[1]) + "-" + str(i[2]))
        real_frag = i[3].split(',')

        lift_real[lifted_frag] = []
        for frag in real_frag:
            lift_real[lifted_frag].append(frag)

real_lift = {}
with open("../../result/alignments/mouse2human/AlignmentStatistics_TBA_mouse2human_withoutnull.txt") as f2:
    for i in f2.readlines()[1:]:
        i = i.strip("\n")
        i = i.split("\t")
        real = i[0].split(":")
        frag = (str(real[0]) + ":" + str(real[1]) + "-" + str(real[2]))

        lift = i[1].split(":")
        lifted_frag = (str(lift[0]) + ":" + str(lift[1]) + "-" + str(lift[2]))

        real_lift[lifted_frag] = frag

# Conservation interaction
output = open("../../result/conservation/mouse2human_conserv_interaction_simul.txt", 'w')
if os.stat("../../result/conservation/mouse2human_conserv_interaction_simul.txt").st_size == 0:
    output.write("bait_mouse\tPIR_mouse\tbait_lift\tPIR_lift\tbait_human\tPIR_human\n")

conserv_interaction = conserv_bait = conserv_PIR = 0
for lifted_bait in inter_lift.keys():
    if lifted_bait in lift_real.keys():
        conserv_bait += 1
        real_bait = lift_real[lifted_bait]

        for lifted_PIR in inter_lift[lifted_bait]:
            if lifted_PIR in lift_real.keys():
                conserv_PIR += 1
                real_PIR = lift_real[lifted_PIR]

                for bait in real_bait:
                    if bait in inter_real.keys():
                        for PIR in real_PIR:
                            if PIR in inter_real[bait]:
                                conserv_interaction += 1
                                output.write(real_lift[lifted_bait]+'\t'+real_lift[lifted_PIR]+'\t' +
                                             lifted_bait+'\t'+lifted_PIR+'\t'+bait+'\t'+PIR+'\n')


print("Conserv baits: ", conserv_bait)
print("Conserv PIRs: ", conserv_PIR)
print("Conserv interactions: ", conserv_interaction)
print("Lifted interactions: ", len(dist_lift))
print("Lifted interactions: ", len(dist_real))

output.close()
