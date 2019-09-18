#!/usr/bin/env python3
# coding=utf-8

import os
import numpy as np


def conserv_seq(origin_sp, target_sp, data):
    print(origin_sp, "to", target_sp, ":", data)
    # Alignment scores
    def dic_score(file):
        frag_conserv = {}
        with open("../../result/alignments/"+file) as f2:
            for i in f2.readlines()[1:]:
                i = i.strip("\n")
                i = i.split("\t")
                frag = i[0].split(':')
                frag = str(frag[0]+':' + str(int(frag[1])+1) + ':' + str(frag[2]))
                all_ungapped = int(i[4])
                all_identical = int(i[5])
                all_length = int(i[8])
                exclude_ungapped = int(i[6])
                exclude_identical = int(i[7])
                all_exclude = int(i[9])

                frag_conserv[frag] = np.array([all_ungapped, all_identical, all_length, exclude_ungapped,
                                               exclude_identical, all_exclude])

        return frag_conserv

    all_exon = "align_stats/AlignmentStatistics_PECAN_ExcludingExons_"+origin_sp+"2"+target_sp+".txt"
    score_extract_all = dic_score(all_exon)
    # all_genes = "AlignmentStatistics_PECAN_ExcludingGenes_"+origin_sp+"2"+target_sp+".txt"
    # score_extract_genes = dic_score(all_genes)
    # coding_exon = "AlignmentStatistics_pecan_Excluding_coding_Exons.txt"
    # score_extract_coding = dic_score(coding_exon)
    # nocoding_exon = "AlignmentStatistics_pecan_Excluding_nocoding_Exons.txt"
    # score_extract_nocoding = dic_score(nocoding_exon)

    print("Alignment infos done ! ")

    # Composition in base pairs
    def dic_pb(file):
        elem_pb = {}
        with open("../../data/"+origin_sp+"/overlap/"+origin_sp+file) as f2:
            for i in f2.readlines()[1:]:
                i = i.strip("\n")
                i = i.split("\t")
                frag = ('chr' + str(i[0]) + ':' + str(i[1]) + ':' + str(i[2]))
                pb = int(i[5])
                elem_pb[frag] = pb

        return elem_pb

    all = "_frag_overlap_all_exons.txt"
    all_exon = dic_pb(all)
    coding = "_frag_overlap_coding_exons.txt"
    coding_exon = dic_pb(coding)
    nocoding = "_frag_overlap_nocoding_exons.txt"
    nocoding_exon = dic_pb(nocoding)
    repeat = "_frag_overlap_repeatmasker.txt"
    repeat_pb = dic_pb(repeat)
    phastcons = "_frag_overlap_phastcons_noexonic250.txt"
    phastcons_pb = dic_pb(phastcons)

    # Composition in count number
    def dic_count(file):
        elem_pb = {}
        with open("../../data/"+origin_sp+"/overlap/"+origin_sp+file) as f2:
            for i in f2.readlines()[1:]:
                i = i.strip("\n")
                i = i.split("\t")
                frag = ('chr' + str(i[0]) + ':' + str(i[1]) + ':' + str(i[2]))

                if i[3] == "NA":
                    count = 0
                else:
                    count = len(i[3].split(","))

                elem_pb[frag] = count

        return elem_pb

    TSS = "_frag_overlap_TSS.txt"
    TSS_count = dic_count(TSS)
    CAGE = "_frag_overlap_CAGE.txt"
    CAGE_count = dic_count(CAGE)
    if origin_sp == "human":
        ENCODE = "_frag_overlap_ENCODE.txt"
        ENCODE_count = dic_count(ENCODE)
        GRO_seq = "_frag_overlap_GRO_seq.txt"
        GRO_seq_count = dic_count(GRO_seq)
        RoadMap = "_frag_overlap_RoadMap.txt"
        RoadMap_count = dic_count(RoadMap)

    print("Calcul overlap done ! ")

    # Duplication score
    frag_dupli = {}
    with open("../../data/"+origin_sp+"/"+origin_sp+"_restriction_fragments_duplication_0.8.txt") as f1:
        for i in f1.readlines()[1:]:
            i = i.strip("\n")
            i = i.split("\t")
            frag = i[0].split(':')
            frag = str(frag[0]+':' + str(int(frag[1])+1) + ':' + str(frag[2]))
            frag_dupli[frag] = int(i[3])-1

    print("Score dupli done ! ")

    # Matching with contacted fragments
    inter = {}
    inter_bait = {}
    contact_unbaited = {}
    PIR_baited = {}
    CAGE_contact = {}
    RoadMap_contact = {}
    ENCODE_contact = {}
    GRO_seq_contact = {}
    link = {}
    link_bait = {}
    if data == "_simul":
        infile = "/Simulations/simulations_"+origin_sp+"_10Mb_bin5kb_fragoverbin_chr"+merge+".txt"
    else:
        infile = "/all_interactions/all_interactions_chr"+merge+".txt"

    not_conserv = not_dupli = 0
    with open("../../data/"+origin_sp+infile) as f3:
        for i in f3.readlines()[1:]:
            i = i.strip("\n")
            i = i.split("\t")
            bait = (str(i[0]) + ":" + str(i[1]) + ":" + str(i[2]))
            merged_PIR = (str(i[3]) + ":" + str(i[4]) + ":" + str(i[5]))

            if i[0] == i[3]:
                dist_obs = float(i[7])
                if 25000 < dist_obs <= 10000000:
                    PIR_baited[merged_PIR] = i[6]

                    ## Bait side
                    if bait in link_bait.keys():
                        link_bait[bait].append(merged_PIR)
                    else:
                        link_bait[bait] = [merged_PIR]
                    if bait in inter_bait.keys():
                        inter_bait[bait].append(dist_obs)
                    else:
                        inter_bait[bait] = [dist_obs]

                    if bait not in score_extract_all.keys():
                        not_conserv += 1
                        score_extract_all[bait] = np.array([0, 0, 0, 0, 0, 0])
                    if bait not in frag_dupli.keys():
                        not_dupli += 1
                        frag_dupli[bait] = 0

                    ## PIR side
                    if merged_PIR in link.keys():
                        link[merged_PIR].append(bait)
                    else:
                        link[merged_PIR] = [bait]
                    if merged_PIR in inter.keys():
                        inter[merged_PIR].append(dist_obs)
                    else:
                        inter[merged_PIR] = [dist_obs]


                    # Adding PIR with no information
                    for PIR in i[8].split(','):
                        if PIR not in score_extract_all.keys():
                            not_conserv += 1
                            score_extract_all[PIR] = np.array([0, 0, 0, 0, 0, 0])
                        if PIR not in frag_dupli.keys():
                            not_dupli += 1
                            frag_dupli[PIR] = 0

                    if len(i[8].split(',')) > 1:
                        for PIR in i[8].split(','):
                        # Merging info of merged fragment
                            if merged_PIR not in score_extract_all.keys():
                                score_extract_all[merged_PIR] = score_extract_all[PIR]
                            else:
                                score_extract_all[merged_PIR] += score_extract_all[PIR]

                        frag_dupli[merged_PIR] = sum(frag_dupli[PIR] for PIR in i[8].split(','))
                        all_exon[merged_PIR] = sum(all_exon[PIR] for PIR in i[8].split(','))
                        coding_exon[merged_PIR] = sum(coding_exon[PIR] for PIR in i[8].split(','))
                        nocoding_exon[merged_PIR] = sum(nocoding_exon[PIR] for PIR in i[8].split(','))
                        repeat_pb[merged_PIR] = sum(repeat_pb[PIR] for PIR in i[8].split(','))
                        phastcons_pb[merged_PIR] = sum(phastcons_pb[PIR] for PIR in i[8].split(','))
                        TSS_count[merged_PIR] = sum(TSS_count[PIR] for PIR in i[8].split(','))
                        CAGE_count[merged_PIR] = sum(CAGE_count[PIR] for PIR in i[8].split(','))

                        if origin_sp == "human":
                            ENCODE_count[merged_PIR] = sum(ENCODE_count[PIR] for PIR in i[8].split(','))
                            GRO_seq_count[merged_PIR] = sum(GRO_seq_count[PIR] for PIR in i[8].split(','))
                            RoadMap_count[merged_PIR] = sum(RoadMap_count[PIR] for PIR in i[8].split(','))

                    if bait not in contact_unbaited.keys():
                        contact_unbaited[bait] = []
                        CAGE_contact[bait] = 0
                        RoadMap_contact[bait] = 0
                        ENCODE_contact[bait] = 0
                        GRO_seq_contact[bait] = 0

                    if i[6] == "unbaited":
                        contact_unbaited[bait].append(merged_PIR)
                    if CAGE_count[merged_PIR] > 0:
                        CAGE_contact[bait] += 1
                    if origin_sp == "human":
                        if RoadMap_count[merged_PIR] > 0:
                            RoadMap_contact[bait] += 1
                        if ENCODE_count[merged_PIR] > 0:
                            ENCODE_contact[bait] += 1
                        if GRO_seq_count[merged_PIR] > 0:
                            GRO_seq_contact[bait] += 1

    print("Not conserved frag:", not_conserv, "; No dupli info:", not_dupli, 'over', len(inter.keys()), 'contacted frag')

    print("Writting output...")
    output = open("../../result/alignments/PIR_cons_all_overlap_PECAN_"+origin_sp+"2"+target_sp+data+merge+".txt", 'w')
    output_bait = open("../../result/alignments/Bait_cons_all_overlap_PECAN_" + origin_sp + "2" + target_sp + data + merge+".txt",'w')

    if origin_sp == "human":
        if os.stat("../../result/alignments/PIR_cons_all_overlap_PECAN_"+origin_sp+"2"+target_sp+data+merge+".txt").st_size == 0:
            output.write("chr\tstart\tend\tCAGE_count\tENCODE_count\tGRO_seq_count\tRoadMap_count\tnb_PIRcontact"
                         "\tnb_baitcontact\tmidist_obs\tbaited\tall_ungapped\tall_identical\tall_length"
                         "\texclude_ungapped\texclude_identical\tall_exclude\tduplication\tall_exon_pb"
                         "\tcoding_exon_pb\tnocoding_exon_pb\trepeat_pb\tphastcons_noexonic250\tTSS_count\n")

        for PIR in inter.keys():
            bait_contact = [len(inter_bait[bait]) for bait in link[PIR]]
            if PIR not in score_extract_all.keys():
                score_extract_all[PIR] = np.array([0, 0, 0, 0, 0, 0])
            if PIR not in frag_dupli.keys():
                frag_dupli[PIR] = 0

            output.write(PIR.split(':')[0] + '\t' + PIR.split(':')[1] + '\t' + PIR.split(':')[2]
                         + '\t' + str(CAGE_count[PIR]) + '\t' + str(ENCODE_count[PIR])
                         + '\t' + str(GRO_seq_count[PIR]) + '\t' + str(RoadMap_count[PIR])
                         + '\t' + str(len(inter[PIR])) + '\t' + str(np.mean(bait_contact)) + '\t' + str(np.median(inter[PIR]))
                         + '\t' + str(PIR_baited[PIR]) + '\t' + str(score_extract_all[PIR][0])
                         + '\t' + str(score_extract_all[PIR][1]) + '\t' + str(score_extract_all[PIR][2])
                         + '\t' + str(score_extract_all[PIR][3]) + '\t' + str(score_extract_all[PIR][4])
                         + '\t' + str(score_extract_all[PIR][5]) + '\t' + str(frag_dupli[PIR])
                         + '\t' + str(all_exon[PIR]) + '\t' + str(coding_exon[PIR]) + '\t' + str(nocoding_exon[PIR])
                         + '\t' + str(repeat_pb[PIR]) + '\t' + str(phastcons_pb[PIR])
                         + '\t' + str(TSS_count[PIR]) + '\n')

        if os.stat("../../result/alignments/Bait_cons_all_overlap_PECAN_"+origin_sp+"2"+target_sp+data+merge+".txt").st_size == 0:
            output_bait.write("chr\tstart\tend\tCAGE_contact\tENCODE_contact\tGRO_seq_contact\tRoadMap_contact"
                              "\tnb_Baitcontact\tnb_unbaited_contact\tnb_PIRcontact\tmidist_obs\tall_ungapped"
                              "\tall_identical\tall_length\texclude_ungapped\texclude_identical\tall_exclude"
                              "\tduplication\tall_exon_pb\tcoding_exon_pb\tnocoding_exon_pb\trepeat_pb"
                              "\tphastcons_noexonic250\tTSS_count\n")

        for Bait in inter_bait.keys():
            PIR_contact = [len(inter[PIR]) for PIR in link_bait[Bait]]

            output_bait.write(Bait.split(':')[0] + '\t' + Bait.split(':')[1] + '\t' + Bait.split(':')[2]
                              + '\t' + str(CAGE_contact[Bait]) + '\t' + str(ENCODE_contact[Bait])
                              + '\t' + str(GRO_seq_contact[Bait]) + '\t' + str(RoadMap_contact[Bait])
                              + '\t' + str(len(inter_bait[Bait])) + '\t' + str(len(contact_unbaited[Bait]))
                              + '\t' + str(np.mean(PIR_contact)) + '\t' + str(np.median(inter_bait[Bait]))
                              + '\t' + str(score_extract_all[Bait][0]) + '\t' + str(score_extract_all[Bait][1])
                              + '\t' + str(score_extract_all[Bait][2]) + '\t' + str(score_extract_all[Bait][3])
                              + '\t' + str(score_extract_all[Bait][4]) + '\t' + str(score_extract_all[Bait][5])
                              + '\t' + str(frag_dupli[Bait]) + '\t' + str(all_exon[Bait])
                              + '\t' + str(coding_exon[Bait]) + '\t' + str(nocoding_exon[Bait])
                              + '\t' + str(repeat_pb[Bait]) + '\t' + str(phastcons_pb[Bait]) + '\t' + str(TSS_count[Bait]) + '\n')
    else:
        if os.stat("../../result/alignments/PIR_cons_all_overlap_PECAN_"+origin_sp+"2"+target_sp+data+merge+".txt").st_size == 0:
            output.write("chr\tstart\tend\tCAGE_count\tnb_PIRcontact\tnb_baitcontact\tmidist_obs\baited\ttall_ungapped\tall_identical\tall_length"
                         "\texclude_ungapped\texclude_identical\tall_exclude\tduplication\tall_exon_pb\tcoding_exon_pb"
                         "\tnocoding_exon_pb\trepeat_pb\tphastcons_noexonic250\tTSS_count\n")

        for PIR in inter.keys():
            bait_contact = [len(inter_bait[bait]) for bait in link[PIR]]
            if PIR not in score_extract_all.keys():
                score_extract_all[PIR] = np.array([0, 0, 0, 0, 0, 0])
            if PIR not in frag_dupli.keys():
                frag_dupli[PIR] = 0
            output.write(PIR.split(':')[0] + '\t' + PIR.split(':')[1] + '\t' + PIR.split(':')[2]
                         + '\t' + str(CAGE_count[PIR]) + '\t' + str(len(inter[PIR]))
                         + '\t' + str(np.mean(bait_contact)) + '\t' + str(np.median(inter[PIR]))
                         + '\t' + str(PIR_baited[PIR])
                         + '\t' + str(score_extract_all[PIR][0]) + '\t' + str(score_extract_all[PIR][1])
                         + '\t' + str(score_extract_all[PIR][2]) + '\t' + str(score_extract_all[PIR][3])
                         + '\t' + str(score_extract_all[PIR][4]) + '\t' + str(score_extract_all[PIR][5])
                         + '\t' + str(frag_dupli[PIR]) + '\t' + str(all_exon[PIR]) + '\t' + str(coding_exon[PIR])
                         + '\t' + str(nocoding_exon[PIR]) + '\t' + str(repeat_pb[PIR]) + '\t' + str(phastcons_pb[PIR])
                         + '\t' + str(TSS_count[PIR]) + '\n')

        if os.stat("../../result/alignments/Bait_cons_all_overlap_PECAN_"+origin_sp+"2"+target_sp+data+merge+".txt").st_size == 0:
            output_bait.write("chr\tstart\tend\tCAGE_count\tnb_contact\ttnb_unbaited_contact\tnb_PIRcontact\tmidist_obs\tall_ungapped\tall_identical\tall_length"
                              "\texclude_ungapped\texclude_identical\tall_exclude\tduplication\tall_exon_pb\tcoding_exon_pb"
                              "\tnocoding_exon_pb\trepeat_pb\tphastcons_noexonic250\tTSS_count\n")

        for Bait in inter_bait.keys():
            PIR_contact = [len(inter[PIR]) for PIR in link_bait[Bait]]
            output_bait.write(Bait.split(':')[0] + '\t' + Bait.split(':')[1] + '\t' + Bait.split(':')[2]
                              + '\t' + str(CAGE_count[Bait]) + '\t' + str(len(inter_bait[Bait]))
                              + '\t' + str(len(contact_unbaited[Bait])) + '\t' + str(np.mean(PIR_contact)) + '\t' +
                              str(np.median(inter_bait[Bait])) + '\t' + str(score_extract_all[Bait][0]) + '\t' +
                              str(score_extract_all[Bait][1]) + '\t' + str(score_extract_all[Bait][2]) + '\t' +
                              str(score_extract_all[Bait][3]) + '\t' + str(score_extract_all[Bait][4]) + '\t' +
                              str(score_extract_all[Bait][5]) + '\t' + str(frag_dupli[Bait]) + '\t' +
                              str(all_exon[Bait]) + '\t' + str(coding_exon[Bait]) + '\t' + str(nocoding_exon[Bait]) + '\t' +
                              str(repeat_pb[Bait]) + '\t' + str(phastcons_pb[Bait]) + '\t' + str(TSS_count[Bait]) + '\n')

    output.close()


origin_sp = "mouse"
target_sps = ["human"] #, "dog", "cow", "elephant", "opossum", "chicken"]
datas = ["", "_simul"]
merge = "_merged"

for sp in target_sps:
    for data in datas:
        conserv_seq(origin_sp, sp, data)

print("All done ! ")

