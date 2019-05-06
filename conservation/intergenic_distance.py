#!/usr/bin/env python3
# coding=utf-8

import os

ref = "mouse"
sp = "human"

ortho = {}
gen_target = {}
gen_ref = {}
gen_sp = {}

with open("../../data/ensembl_ortho/Ortho_"+ref+"_"+sp+"_Ensembl91_one2one.txt", 'r') as f1:
    for i in f1.readlines():
        i = i.split('\t')
        ortho[i[0]] = i[4]

        chr_ref = str(i[1])
        mid_gene_ref = (int(i[3])+int(i[2]))/2
        if chr_ref not in gen_ref.keys():
            gen_ref[chr_ref] = [(i[0], mid_gene_ref)]
        else:
            gen_ref[chr_ref].append((i[0], mid_gene_ref))

        chr_sp = str(i[5])
        mid_gene_sp = (int(i[7]) + int(i[6])) / 2
        gen_target[i[4]] = (chr_sp, mid_gene_sp)

print('Nombre de gène orthologue:', sum(len(gen_ref[i]) for i in gen_ref.keys()))


output = open("../../result/conservation/intergenic_distance_" + ref + "2" + sp + ".txt", 'w')
if os.stat("../../result/conservation/intergenic_distance_" + ref + "2" + sp + ".txt").st_size == 0:
    output.write("gene_pair\torigin_dist\ttarget_dist\n")

paires = bad_paires = paires_ok = 0
inter_dist = {}

for chr in gen_ref.keys():
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
                        output.write(str(gene[0]) + "_" + str(next_gen[0]) + '\t' +
                                     str(dist_ref) + "\t" + str(dist_target) + "\n")
                        paires_ok += 1

                else:
                    bad_paires += 1

                index_next_gen += 1

output.close()

print('Nombre total de paires de gènes valides:', paires)
print('Nombre de paires de gènes ok:', paires_ok)
print('Nombre de paires de gènes non valides:', bad_paires)
