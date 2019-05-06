#!/usr/bin/env python3
# coding=utf-8

import os

ref = "mouse"
sp = "human"
path = "/home/laverre/Documents/Regulatory_Landscape/"

ortho = {}
gen_target = {}
gen_ref = {}
gen_sp = {}

with open(path+"result/alignments/"+ref+"2"+sp+"/AlignmentStatistics_TBA_"+ref+"2"+sp+"_withoutnull.txt", 'r') as f1:
    for i in f1.readlines()[1:]:
        i = i.split('\t')
        ortho[i[0]] = i[1]
        frag_ref = i[0].split(':')

        chr_ref = str(frag_ref[0])
        mid_gene_ref = (int(frag_ref[1])+int(frag_ref[2]))/2
        if chr_ref not in gen_ref.keys():
            gen_ref[chr_ref] = [(i[0], mid_gene_ref)]
        else:
            gen_ref[chr_ref].append((i[0], mid_gene_ref))

        frag_sp = i[1].split(':')
        chr_sp = str(frag_sp[0])
        mid_gene_sp = (int(frag_sp[1]) + int(frag_sp[2])) / 2
        gen_target[i[1]] = (chr_sp, mid_gene_sp)

print('Nombre de gène orthologue:', sum(len(gen_ref[i]) for i in gen_ref.keys()))


output = open(path+"result/conservation/interfrag_distance_" + ref + "2" + sp + ".txt", 'w')
if os.stat(path+"result/conservation/interfrag_distance_" + ref + "2" + sp + ".txt").st_size == 0:
    output.write("origin_dist\ttarget_dist\n")

paires = bad_paires = paires_ok = 0
inter_dist = {}

count_chr = 0
for chr in gen_ref.keys():
    count_chr += 1
    print("Chromosome :", count_chr, "sur", len(gen_ref.keys()))
    if len(gen_ref[chr]) > 1:  # chromo contient + de 1 gène
        for gene in gen_ref[chr]:
            index_next_gen = gen_ref[chr].index(gene)+1  # première paire de gène du chr

            while index_next_gen < len(gen_ref[chr]):
                next_gen = gen_ref[chr][index_next_gen]
                gen_ortho = gen_target[ortho[gene[0]]]
                gen_next_ortho = gen_target[ortho[next_gen[0]]]

                if gen_ortho[0] == gen_next_ortho[0]:  # paires homologues sur le même chromo
                    dist_ref = abs(next_gen[1] - gene[1])
                    dist_target = abs(gen_ortho[1] - gen_next_ortho[1])
                    paires += 1

                    if 25000 < dist_ref < 10000000:
                        output.write(str(dist_ref) + "\t" + str(dist_target) + "\n")
                        paires_ok += 1

                else:
                    bad_paires += 1

                index_next_gen += 1

output.close()

print('Nombre total de paires de gènes valides:', paires)
print('Nombre de paires de gènes ok:', paires_ok)
print('Nombre de paires de gènes non valides:', bad_paires)
