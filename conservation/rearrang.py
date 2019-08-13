#!/usr/bin/env python3
# coding=utf-8

origin_sp = "human"
#target_sps = ["human", "dog", "cow", "elephant", "opossum", "chicken"]
#datas = ["_simul", ""]  # or "_simul"
#+ origin_sp + infile

conserv = {}
with open("../../result/syntenie/human2mouse_converted_fragments_duplication.txt") as f3:
    for i in f3.readlines()[1:]:
        i = i.strip("\n")
        i = i.split("\t")
        frag = i[0]
        coord = frag.split(':')
        conserv[coord[0]+':'+coord[1]+'-'+coord[2]] = 0

print('Conserv list ok!')

bloc = {}
with open("../../result/syntenie/blocs_gros") as f3:
    for i in f3.readlines()[1:]:
        i = i.strip("\n")
        i = i.split(" ")
        bloc[i[0]] = i[1]

print('Bloc dic ok!')

same_bloc = bloc_diff = total = bait_ok = PIR_ok = total_conserv = 0
PIR_list = {}
bait_list = {}
infile = "/all_interactions/all_interactions_chr.txt"
#infile =  "/Simulations/simulations_"+origin_sp+"_10Mb_bin5kb_fragoverbin_chr.txt"
with open("../../data/" + origin_sp + infile) as f3:
    for i in f3.readlines()[1:]:
        i = i.strip("\n")
        i = i.split("\t")
        bait = (str(i[0]) + ":" + str(int(i[1])-1) + "-" + str(i[2]))
        PIR = (str(i[3]) + ":" + str(int(i[4])-1) + "-" + str(i[5]))
        midbait = ((int(i[1]) + int(i[2])) / 2)
        midcontact = ((int(i[4]) + int(i[5])) / 2)

        if i[0] == i[3]:
            if 25000 < abs(midbait - midcontact) <= 10000000:
                total += 1
                bait_list[bait] = 0
                PIR_list[PIR] = 0

                if bait in conserv.keys():
                    if PIR in conserv.keys():
                        total_conserv += 1

                        if bait in bloc.keys():
                            if PIR in bloc.keys():
                                if bloc[bait] == bloc[PIR]:
                                    same_bloc += 1
                                else:
                                    bloc_diff += 1

print("Total interactions :", total)
print("Total interactions conservées :", total_conserv)
print("Nombre interactions même bloc:", same_bloc)
print("Nombre interactions bloc diff:", bloc_diff)
print("Nombre total de bait:", len(bait_list.keys()))
print("Nombre total de PIR:", len(PIR_list.keys()))
