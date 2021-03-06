#!/usr/bin/env python3
# coding=utf-8

import os
import numpy as np


def conserv_seq(origin_sp, target_sp, data):
    print(origin_sp, "to", target_sp, ":", data)
    # Alignment scores
    def dic_score(file):
        frag_conserv = {}
        with open("../../result/alignments/align_stats/"+file) as f2:
            for i in f2.readlines()[1:]:
                i = i.strip("\n")
                i = i.split("\t")
                frag = i[0].split(':')
                frag = str(frag[0]+':' + str(int(frag[1])+1) + ':' + str(frag[2]))

                all_ungapped = int(i[4]) / (int(i[8])+1)
                all_identical = int(i[5]) / (int(i[8])+1)

                exclude_ungapped = int(i[6]) / (int(i[9])+1)
                exclude_identical = int(i[7]) / (int(i[9])+1)

                if file == all_exon:
                    frag_conserv[frag] = (str(all_ungapped) + '\t' + str(all_identical) + '\t' +
                                          str(exclude_ungapped) + '\t' + str(exclude_identical))
                else:
                    frag_conserv[frag] = (str(exclude_ungapped))

        return frag_conserv

    all_exon = "AlignmentStatistics_PECAN_ExcludingExons_"+origin_sp+"2"+target_sp+".txt"
    score_extract_all = dic_score(all_exon)
    all_genes = "AlignmentStatistics_PECAN_ExcludingGenes_"+origin_sp+"2"+target_sp+".txt"
    score_extract_genes = dic_score(all_genes)
    # coding_exon = "AlignmentStatistics_pecan_Excluding_coding_Exons.txt"
    # score_extract_coding = dic_score(coding_exon)
    # nocoding_exon = "AlignmentStatistics_pecan_Excluding_nocoding_Exons.txt"
    # score_extract_nocoding = dic_score(nocoding_exon)

    print("Calcul alignment scores done ! ")

    # Composition in base pairs
    def dic_pb(file):
        elem_pb = {}
        with open("../../data/"+origin_sp+"/overlap/"+origin_sp+file) as f2:
            for i in f2.readlines()[1:]:
                i = i.strip("\n")
                i = i.split("\t")
                frag = ('chr' + str(i[0]) + ':' + str(i[1]) + ':' + str(i[2]))

                length = str(int(i[4]))
                pb = str(int(i[5]))

                if file == "_frag_overlap_all_exons.txt":
                    elem_pb[frag] = str(length + '\t' + pb)
                else:
                    elem_pb[frag] = str(pb)

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

                elem_pb[frag] = str(count)

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
            frag_dupli[frag] = str(int(i[3])-1)

    print("Score dupli done ! ")

    # Matching with contacted fragments
    inter = {}
    inter_bait = {}
    link = {}
    link_bait = {}
    if data == "_simul":
        infile = "/Simulations/simulations_"+origin_sp+"_10Mb_bin5kb_fragoverbin_chr.txt"
    else:
        infile = "/all_interactions/all_interactions_chr.txt"

    not_conserv = not_dupli = 0
    with open("../../data/"+origin_sp+infile) as f3:
        for i in f3.readlines()[1:]:
            i = i.strip("\n")
            i = i.split("\t")
            midbait = ((int(i[1]) + int(i[2])) / 2)
            midcontact = ((int(i[4]) + int(i[5])) / 2)
            bait = (str(i[0]) + ":" + str(i[1]) + ":" + str(i[2]))
            PIR = (str(i[3]) + ":" + str(i[4]) + ":" + str(i[5]))
            dist_obs = abs(midbait - midcontact)

            if i[0] == i[3]:
                if 25000 < abs(midbait - midcontact) <= 10000000:

                    ### PIR side
                    if PIR in link.keys():
                        link[PIR].append(bait)
                    else:
                        link[PIR] = [bait]
                    if PIR in inter.keys():
                        inter[PIR].append(dist_obs)
                    else:
                        inter[PIR] = [dist_obs]
                    if PIR not in score_extract_all.keys():
                        not_conserv += 1
                        score_extract_all[PIR] = (str(0) + '\t' + str(0) + '\t' + str(0) + '\t' + str(0))
                    if PIR not in score_extract_genes.keys():
                        score_extract_genes[PIR] = str(0)
                    if PIR not in frag_dupli.keys():
                        not_dupli += 1
                        frag_dupli[PIR] = 'NA'

                    ## Bait side
                    if bait in link_bait.keys():
                        link_bait[bait].append(PIR)
                    else:
                        link_bait[bait] = [PIR]
                    if bait in inter_bait.keys():
                        inter_bait[bait].append(dist_obs)
                    else:
                        inter_bait[bait] = [dist_obs]
                    if bait not in score_extract_all.keys():
                        not_conserv += 1
                        score_extract_all[bait] = (str(0) + '\t' + str(0) + '\t' + str(0) + '\t' + str(0))
                    if bait not in score_extract_genes.keys():
                        score_extract_genes[bait] = str(0)
                    if bait not in frag_dupli.keys():
                        not_dupli += 1
                        frag_dupli[bait] = 'NA'

    print("Not conserved frag:", not_conserv, "; No dupli info:", not_dupli, 'over', len(inter.keys()), 'contacted frag')


    print("Writting output...")
    output = open("../../result/alignments/PIR_cons_all_overlap_PECAN_"+origin_sp+"2"+target_sp+data+".txt3", 'w')
    output_bait = open("../../result/alignments/Bait_cons_all_overlap_PECAN_" + origin_sp + "2" + target_sp + data + ".txt3",'w')

    if origin_sp == "human":
        if os.stat("../../result/alignments/PIR_cons_all_overlap_PECAN_"+origin_sp+"2"+target_sp+data+".txt3").st_size == 0:
            output.write("chr\tstart\tend\tCAGE_count\tENCODE_count\tGRO_seq_count\tRoadMap_count\tnb_PIRcontact"
                         "\tnb_baitcontact\tmidist_obs\tTotal_ungapped\tTotal_identical\tAllexon_ungapped\tAllexon_identical"
                         "\tAllgene_ungapped\tduplication\tlength\tall_exon_pb\tcoding_exon_pb\tnocoding_exon_pb\trepeat_pb"
                         "\tphastcons_noexonic250\tTSS_count\n")

        for PIR in inter.keys():
            bait_contact = [len(inter_bait[bait]) for bait in link[PIR]]
            output.write(PIR.split(':')[0] + '\t' + PIR.split(':')[1] + '\t' + PIR.split(':')[2]
                         + '\t' + CAGE_count[PIR] + '\t' + ENCODE_count[PIR] + '\t' + GRO_seq_count[PIR]
                         + '\t' + RoadMap_count[PIR] + '\t' + str(len(inter[PIR])) + '\t' + str(np.mean(bait_contact))
                         + '\t' + str(np.median(inter[PIR])) + '\t' + score_extract_all[PIR]
                         + '\t' + score_extract_genes[PIR] + '\t' + frag_dupli[PIR] + '\t' + all_exon[PIR]
                         + '\t' + coding_exon[PIR] + '\t' + nocoding_exon[PIR] + '\t' + repeat_pb[PIR]
                         + '\t' + phastcons_pb[PIR] + '\t' + TSS_count[PIR] + '\n')

        if os.stat("../../result/alignments/Bait_cons_all_overlap_PECAN_"+origin_sp+"2"+target_sp+data+".txt3").st_size == 0:
            output_bait.write("chr\tstart\tend\tCAGE_count\tENCODE_count\tGRO_seq_count\tRoadMap_count\tnb_Baitcontact"
                              "\tnb_PIRcontact\tmidist_obs\tTotal_ungapped\tTotal_identical\tAllexon_ungapped"
                              "\tAllexon_identical\tAllgene_ungapped\tduplication\tlength\tall_exon_pb\tcoding_exon_pb"
                              "\tnocoding_exon_pb\trepeat_pb\tphastcons_noexonic250\tTSS_count\n")

        for Bait in inter_bait.keys():
            PIR_contact = [len(inter[PIR]) for PIR in link_bait[Bait]]
            output_bait.write(Bait.split(':')[0] + '\t' + Bait.split(':')[1] + '\t' + Bait.split(':')[2]
                              + '\t' + CAGE_count[Bait] + '\t' + ENCODE_count[Bait] + '\t' + GRO_seq_count[Bait]
                              + '\t' + RoadMap_count[Bait] + '\t' + str(len(inter_bait[Bait])) + '\t' + str(np.mean(PIR_contact))
                              + '\t' + str(np.median(inter_bait[Bait])) + '\t' + score_extract_all[Bait]
                              + '\t' + score_extract_genes[Bait] + '\t' + frag_dupli[Bait] + '\t' + all_exon[Bait]
                              + '\t' + coding_exon[Bait] + '\t' + nocoding_exon[Bait] + '\t' + repeat_pb[Bait]
                              + '\t' + phastcons_pb[Bait] + '\t' + TSS_count[Bait] + '\n')
    else:
        if os.stat("../../result/alignments/PIR_cons_all_overlap_PECAN_"+origin_sp+"2"+target_sp+data+".txt3").st_size == 0:
            output.write("chr\tstart\tend\tCAGE_count\tnb_contact\tmidist_obs\tTotal_ungapped\tTotal_identical\tAllexon_ungapped\tAllexon_identical"
                         "\tAllgene_ungapped\tduplication\tlength\tall_exon_pb\tcoding_exon_pb\tnocoding_exon_pb\trepeat_pb"
                         "\tphastcons_noexonic250\tTSS_count\n")

        for PIR in inter.keys():
            output.write(PIR.split(':')[0] + '\t' + PIR.split(':')[1] + '\t' + PIR.split(':')[2]
                         + '\t' + CAGE_count[PIR]
                         + '\t' + str(len(inter[PIR])) + '\t' + str(np.median(inter[PIR])) + '\t' + score_extract_all[PIR]
                         + '\t' + score_extract_genes[PIR] + '\t' + frag_dupli[PIR] + '\t' + all_exon[PIR] + '\t' + coding_exon[PIR]
                         + '\t' + nocoding_exon[PIR] + '\t' + repeat_pb[PIR] + '\t' + phastcons_pb[PIR] + '\t' + TSS_count[PIR]
                         + '\n')

        if os.stat("../../result/alignments/Bait_cons_all_overlap_PECAN_"+origin_sp+"2"+target_sp+data+".txt3").st_size == 0:
            output_bait.write("chr\tstart\tend\tCAGE_count\tnb_Baitcontact\tnb_PIRcontact\tmidist_obs\tTotal_ungapped"
                              "\tTotal_identical\tAllexon_ungapped\tAllexon_identical\tAllgene_ungapped\tduplication"
                              "\tlength\tall_exon_pb\tcoding_exon_pb\tnocoding_exon_pb\trepeat_pb\tphastcons_noexonic250"
                              "\tTSS_count\n")

        for Bait in inter_bait.keys():
            PIR_contact = [len(inter[PIR]) for PIR in link_bait[Bait]]
            output_bait.write(Bait.split(':')[0] + '\t' + Bait.split(':')[1] + '\t' + Bait.split(':')[2]
                              + '\t' + CAGE_count[Bait]
                              + '\t' + str(len(inter_bait[Bait])) + '\t' + str(np.mean(PIR_contact))
                              + '\t' + str(np.median(inter_bait[Bait])) + '\t' + score_extract_all[Bait]
                              + '\t' + score_extract_genes[Bait] + '\t' + frag_dupli[Bait] + '\t' + all_exon[Bait]
                              + '\t' + coding_exon[Bait] + '\t' + nocoding_exon[Bait] + '\t' + repeat_pb[Bait]
                              + '\t' + phastcons_pb[Bait] + '\t' + TSS_count[Bait] + '\n')



    output.close()


origin_sp = "human"
target_sps = ["mouse"] #, "dog", "cow", "elephant", "opossum", "chicken"]
datas = ["_simul", ""]

for sp in target_sps:
    for data in datas:
        conserv_seq(origin_sp, sp, data)

print("All done ! ")