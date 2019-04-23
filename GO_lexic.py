#!/usr/bin/env python3
# coding=utf-8

GO_lexic = {}
with open("../../data/mouse/annotations/GO_lexic.txt", 'r') as f1:
    for i in f1.readlines()[1:]:
        i = i.strip('\n')
        i = i.split('\t')
        GO_lexic[i[0]] = i[1]

Uniprot_lexic = {}
with open("../../data/mouse/annotations/Uniprot_lexic.txt", 'r') as f1:
    for i in f1.readlines()[1:]:
        i = i.strip('\n')
        i = i.split('\t')
        Uniprot_lexic[i[0]] = i[1]

output = open("../../data/mouse/annotations/Gene_developmental_process_annotated_QuickGO.txt", 'w')
with open("../../data/mouse/annotations/Gene_developmental_process", 'r') as f1:
    for line in f1.readlines():
        line = line.strip('\n')
        i = line.split('\t')

        if i[1] in GO_lexic.keys():
            GO = GO_lexic[i[1]]
        else:
            GO = 'NA'

        if i[0] in Uniprot_lexic.keys():
            ID = Uniprot_lexic[i[0]]
            output.write(line + '\t' + ID + '\t' + GO + '\n')








