#!/usr/bin/env python3
# coding=utf-8

import os
import numpy as np

def conserv_seq(origin_sp, target_sp, data):
    print(origin_sp, "to", target_sp, ":", data)

    # Alignment scores
    def dic_score(file):
        frag_conserv = {}
        with open("../../result/alignments/Alignments_statistics/merged/"+file) as f2:
            for i in f2.readlines()[1:]:
                i = i.strip("\n")
                i = i.split("\t")
                frag = i[0].split(':')
                frag = str(frag[0]+':' + str(frag[1]) + ':' + str(frag[2]))
                all_ungapped = int(i[4])
                all_identical = int(i[5])
                all_length = int(i[8])
                exclude_ungapped = int(i[6])
                exclude_identical = int(i[7])
                all_exclude = int(i[9])

                frag_conserv[frag] = np.array([all_ungapped, all_identical, all_length, exclude_ungapped,
                                               exclude_identical, all_exclude])

        return frag_conserv

    all_exon = "AlignmentStatistics_Excluding_all_Exons_"+origin_sp+"2"+target_sp+".txt"
    score_extract_all = dic_score(all_exon)

    # Composition in base pairs
    def dic_pb(file):
        elem_pb = {}
        with open("../../data/"+origin_sp+"/overlap/"+origin_sp+file) as f2:
            for i in f2.readlines()[1:]:
                i = i.strip("\n")
                i = i.split("\t")
                frag = (str(i[1]) + ':' + str(i[2]) + ':' + str(i[3]))
                elem_pb[frag] = 0 if i[4] == "NA" else int(i[6])

        return elem_pb

    all = "_merged_overlap_all_exons_bp.txt"
    all_exon = dic_pb(all)
    all_250 = "_merged_overlap_all_exons250_bp.txt"
    all_exon250 = dic_pb(all_250)
    coding = "_merged_overlap_coding_exons_bp.txt"
    coding_exon = dic_pb(coding)
    nocoding = "_merged_overlap_nocoding_exons_bp.txt"
    nocoding_exon = dic_pb(nocoding)
    repeat = "_merged_overlap_repeatmasker_bp.txt"
    repeat_pb = dic_pb(repeat)
    phastcons = "_merged_overlap_phastcons250_bp.txt"
    phastcons_pb = dic_pb(phastcons)

    # Composition in count number
    def dic_count(file):
        elem_pb = {}
        with open("../../data/"+origin_sp+"/overlap/"+origin_sp+file) as f2:
            for i in f2.readlines()[1:]:
                i = i.strip("\n")
                i = i.split("\t")
                frag = (str(i[1]) + ':' + str(i[2]) + ':' + str(i[3]))
                elem_pb[frag] = 0 if i[4] == "NA" else len(i[4].split(","))

        return elem_pb

    TSS = "_merged_overlap_genes.txt"
    TSS_count = dic_count(TSS)
    CAGE = "_merged_overlap_CAGE.txt"
    CAGE_count = dic_count(CAGE)
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
    link = {}
    link_bait = {}
    contact_unbaited = {}
    PIR_baited = {}
    PIR_nbcell = {}
    CAGE_contact = {}
    RoadMap_contact = {}
    ENCODE_contact = {}
    GRO_seq_contact = {}

    if data == "_simul":
        infile = "/Simulations/simulations_"+origin_sp+"_10Mb_bin5kb_fragoverbin_chr"+merge+".txt"
    else:
        infile = "/all_interactions/all_interactions_chr"+merge+".txt_union_cell"

    with open("../../data/"+origin_sp+infile) as f3:
        for i in f3.readlines()[1:]:
            i = i.strip("\n")
            i = i.split("\t")
            bait = (str(i[0]) + ":" + str(i[1]) + ":" + str(i[2]))
            merged_PIR = (str(i[3]) + ":" + str(i[4]) + ":" + str(i[5]))

            if i[0] == i[3]:
                dist_obs = float(i[7])
                if 25000 <= dist_obs <= 10000000:
                    PIR_baited[merged_PIR] = i[6]
                    PIR_nbcell[merged_PIR] = i[9]

                    # Bait side : names and distances of contacted fragments
                    if bait not in link_bait.keys():
                        link_bait[bait] = [merged_PIR]
                    else:
                        link_bait[bait].append(merged_PIR)
                    if bait not in inter_bait.keys():
                        inter_bait[bait] = [dist_obs]
                    else:
                        inter_bait[bait].append(dist_obs)

                    # PIR side : names and distances of target baits
                    if merged_PIR not in link.keys():
                        link[merged_PIR] = [bait]
                    else:
                        link[merged_PIR].append(bait)
                    if merged_PIR not in inter.keys():
                        inter[merged_PIR] = [dist_obs]
                    else:
                        inter[merged_PIR].append(dist_obs)

                    # Adding non-conserved fragments
                    if bait not in score_extract_all.keys():
                        score_extract_all[bait] = np.array([0, 0, 0, 0, 0, 0])
                    if bait not in frag_dupli.keys():
                        frag_dupli[bait] = 0
                    if merged_PIR not in score_extract_all.keys():
                        score_extract_all[merged_PIR] = np.array([0, 0, 0, 0, 0, 0])
                    if merged_PIR not in frag_dupli.keys():
                        frag_dupli[merged_PIR] = 0

                    # Count contacted enhancer & unbaited contact per bait
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


    print("Writting output...")
    output = open("../../result/conservation/Sequence_conservation/PIR_cons_all_overlap_PECAN_"+origin_sp+"2"+target_sp+data+merge+".txt", 'w')
    #output_bait = open("../../result/alignments/Bait_cons_all_overlap_PECAN_" + origin_sp + "2" + target_sp + data + merge+".txt",'w')

    if os.stat("../../result/conservation/Sequence_conservation/PIR_cons_all_overlap_PECAN_"+origin_sp+"2"+target_sp+data+merge+".txt").st_size == 0:
        output.write("chr\tstart\tend\tCAGE_count\tnb_PIRcontact\tnb_baitcontact\tmidist_obs\tbaited\tnb_cell\ttall_ungapped\tall_identical\tall_length"
                     "\texclude_ungapped\texclude_identical\tall_exclude\tduplication\tall_exon_pb\tall_exon250\tcoding_exon_pb"
                     "\tnocoding_exon_pb\trepeat_pb\tphastcons_noexonic250\tTSS_count\n")

    for PIR in inter.keys():
        bait_contact = [len(inter_bait[bait]) for bait in link[PIR]]

        output.write(PIR.split(':')[0] + '\t' + PIR.split(':')[1] + '\t' + PIR.split(':')[2]
                     + '\t' + str(CAGE_count[PIR]) + '\t' + str(len(inter[PIR]))
                     + '\t' + str(np.mean(bait_contact)) + '\t' + str(np.median(inter[PIR]))
                     + '\t' + str(PIR_baited[PIR]) + '\t' + str(PIR_nbcell[PIR])
                     + '\t' + str(score_extract_all[PIR][0]) + '\t' + str(score_extract_all[PIR][1])
                     + '\t' + str(score_extract_all[PIR][2]) + '\t' + str(score_extract_all[PIR][3])
                     + '\t' + str(score_extract_all[PIR][4]) + '\t' + str(score_extract_all[PIR][5])
                     + '\t' + str(frag_dupli[PIR]) + '\t' + str(all_exon[PIR]) + '\t' + str(all_exon250)
                     + '\t' + str(coding_exon[PIR]) + '\t' + str(nocoding_exon[PIR]) + '\t' + str(repeat_pb[PIR])
                     + '\t' + str(phastcons_pb[PIR]) + '\t' + str(TSS_count[PIR]) + '\n')

    #if os.stat("../../result/alignments/Bait_cons_all_overlap_PECAN_"+origin_sp+"2"+target_sp+data+merge+".txt").st_size == 0:
    #    output_bait.write("chr\tstart\tend\tCAGE_count\tnb_contact\tnb_unbaited_contact\tnb_PIRcontact\tmidist_obs\tall_ungapped\tall_identical\tall_length"
    #                      "\texclude_ungapped\texclude_identical\tall_exclude\tduplication\tall_exon_pb\tcoding_exon_pb"
    #                      "\tnocoding_exon_pb\trepeat_pb\tphastcons_noexonic250\tTSS_count\n")

    #for Bait in inter_bait.keys():
    #    PIR_contact = [len(inter[PIR]) for PIR in link_bait[Bait]]
    #    output_bait.write(Bait.split(':')[0] + '\t' + Bait.split(':')[1] + '\t' + Bait.split(':')[2]
    #                      + '\t' + str(CAGE_count[Bait]) + '\t' + str(len(inter_bait[Bait]))
    #                      + '\t' + str(len(contact_unbaited[Bait])) + '\t' + str(np.mean(PIR_contact)) + '\t' +
    #                      str(np.median(inter_bait[Bait])) + '\t' + str(score_extract_all[Bait][0]) + '\t' +
    #                      str(score_extract_all[Bait][1]) + '\t' + str(score_extract_all[Bait][2]) + '\t' +
    #                      str(score_extract_all[Bait][3]) + '\t' + str(score_extract_all[Bait][4]) + '\t' +
    #                      str(score_extract_all[Bait][5]) + '\t' + str(frag_dupli[Bait]) + '\t' +
    #                      str(all_exon[Bait]) + '\t' + str(coding_exon[Bait]) + '\t' + str(nocoding_exon[Bait]) + '\t' +
    #                      str(repeat_pb[Bait]) + '\t' + str(phastcons_pb[Bait]) + '\t' + str(TSS_count[Bait]) + '\n')

    output.close()


origin_sp = "mouse"
target_sps = ["chicken", "opossum", "rat", "human", "rabbit", "dog", "elephant", "cow", "macaque"]
datas = ["", "_simul"]
merge = "_merged"

for sp in target_sps:
    for data in datas:
        conserv_seq(origin_sp, sp, data)

print("All done ! ")
