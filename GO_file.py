#!/usr/bin/env python3
# coding=utf-8

sp = "human"

GO_lexic = {}
with open("../../data/"+sp+"/annotations/GO_lexic.txt", 'r') as f1:
    for i in f1.readlines()[1:]:
        i = i.strip('\n')
        i = i.split('\t')
        GO_lexic[i[0]] = i[1]

Uniprot_lexic = {}
with open("../../data/"+sp+"/annotations/Uniprot_lexic_immun_process.txt", 'r') as f1:
    for i in f1.readlines()[1:]:
        i = i.strip('\n')
        i = i.split('\t')
        Uniprot_lexic[i[0]] = i[1]

output = open("../../data/"+sp+"/annotations/Gene_immune_system_process_annotated_QuickGO.txt", 'w')
with open("../../data/"+sp+"/annotations/QuickGO_annot_"+sp+"_immun_process.txt", 'r') as f1:
    for i in f1.readlines():
        i = i.strip('\n')
        i = i.split('\t')
        gene_GO = i[3]
        uni_ID = i[1]

        if gene_GO in GO_lexic.keys():
            GO = GO_lexic[gene_GO]
        else:
            GO = 'NA'

        if uni_ID in Uniprot_lexic.keys():
            ID = Uniprot_lexic[uni_ID]
            output.write(uni_ID + '\t' + gene_GO + '\t' + ID + '\t' + GO + '\n')

output.close()

annotation = {}
with open("../../data/"+sp+"/annotations/Gene_immune_system_process_annotated_QuickGO.txt", 'r') as f1:
    for i in f1.readlines():
        i = i.strip('\n')
        i = i.split('\t')
        MGI_ID = i[2]
        annot = i[3]
        if MGI_ID in annotation.keys():
            if annot not in annotation[MGI_ID]:
                annotation[MGI_ID].append(annot)

        else:
            annotation[MGI_ID] = [annot]

output = open("../../data/"+sp+"/annotations/Gene_immune_system_process_annotated_QuickGO_uniq.txt", 'w')
for MGI_ID in annotation.keys():
    count = 0
    output.write(MGI_ID + '\t')

    for annot in annotation[MGI_ID]:
        count += 1
        if count == len(annotation[MGI_ID]):
            output.write(annot + '\n')
        else:
            output.write(annot + ',')

output.close()







