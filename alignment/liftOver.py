#!/usr/bin/env python3
# coding=utf-8
import os

# Séparation bait/enhancer avant liftOver

bait_old = open("../data/promoter_enhancer_contacts/Javierre2016/human/PCHiC_vs_rCHiC_peak_matrix_hg19_bait.txt", "a")
enh_old = open("../data/promoter_enhancer_contacts/Javierre2016/human/PCHiC_vs_rCHiC_peak_matrix_hg19_enh.txt", "a")
dinter = {}
denh = {}
with open("../data/promoter_enhancer_contacts/Javierre2016/human/PCHiC_vs_rCHiC_peak_matrix_hg19.tsv", 'r') as f1:
    for i in f1.readlines()[1:]:
        i = i.split("\t")
        ID_bait = (i[0]+":"+i[1]+":"+i[2])
        ID_enh = (i[5]+":"+i[6]+":"+i[7])
        if ID_bait not in dinter.keys():
            dinter[ID_bait] = [ID_enh]
            bait_old.write(i[0] + "\t" + i[1] + "\t" + i[2] + "\t" + ID_bait + "\n")
        #else:
            #dinter[ID_bait].append(ID_enh)

        if ID_enh not in denh.keys():
            denh[ID_enh] = 1
            enh_old.write(i[5] + "\t" + i[6] + "\t" + i[7] + "\t" + ID_enh + "\n")
"""

# Reformation des interactions en mm10

bait = open("../data/promoter_enhancer_contacts/Schoenfelder2015/mouse/FLC_prom_other_liftOver/FLC_prom_other_bait_mm10.txt", "r")
enh = open("../data/promoter_enhancer_contacts/Schoenfelder2015/mouse/FLC_prom_other_liftOver/FLC_prom_other_enhancer_mm10.txt", "r")
dic_enh = {}
chromo_diff_enh = 0
for i in enh.readlines():
    i = i.strip("\n")
    i = i.split("\t")
    new = (i[0]+":"+i[1]+":"+i[2])
    old_chromo = i[3].split(":")
    if old_chromo[0] == i[0]:
        dic_enh[i[3]] = new
    else :
        chromo_diff_enh +=1

chromo_diff_bait = 0
dic_bait = {}
for i in bait.readlines():
    i = i.strip("\n")
    i = i.split("\t")
    new = (i[0]+":"+i[1]+":"+i[2])
    old = i[3].split(":")
    if old[0] == i[0]:
        dic_bait[i[3]] = new
    else:
        chromo_diff_bait += 1

print(len(dic_enh), len(dic_bait))
print(chromo_diff_enh, chromo_diff_bait)


mm10 = open("../data/promoter_enhancer_contacts/Schoenfelder2015/mouse/FLC_promoter_other_significant_interactions_mm10.txt", 'a')
if os.stat("../data/promoter_enhancer_contacts/Schoenfelder2015/mouse/FLC_promoter_other_significant_interactions_mm10.txt").st_size == 0:
    mm10.write("chr_bait\tstart_bait\tend_bait\texpression_quartile\tchr\tstart\tend\traw_count\tlog(observed/expected)\n")


already = {}
tot = 0
with gzip.open("../data/promoter_enhancer_contacts/Schoenfelder2015/mouse/FLC_promoter_other_significant_interactions_mm9.txt.gz",
    'r') as f1:
    for i in f1.readlines()[1:]:
        tot += 1
        print(tot)
        i = i.decode('UTF-8')
        i = i.split("\t")
        exp = i[5]
        raw = i[9]
        log = i[10]
        old_bait = (i[0]+":"+i[1]+":"+i[2])
        if old_bait in dic_bait.keys():
            bait10 = dic_bait[old_bait]
            bait10 = bait10.split(":")
            chr_bait = bait10[0]
            start_bait = bait10[1]
            end_bait = bait10[2]
            if old_bait not in already.keys():
                already[old_bait] = 1
                for x in dinter[old_bait]:
                    if x in dic_enh.keys():
                        enh10 = dic_enh[x]
                        enh10 = enh10.split(":")
                        chr_enh = enh10[0]
                        start_enh = enh10[1]
                        end_enh = enh10[2]
                        mm10.write(chr_bait+"\t"+start_bait+"\t"+end_bait+"\t"+exp+"\t"+chr_enh+"\t"+start_enh+"\t"
                               +end_enh+"\t"+raw+"\t"+log)
"""

