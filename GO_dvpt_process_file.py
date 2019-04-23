#!/usr/bin/env python3
# coding=utf-8


annotation = {}
with open("../../data/mouse/annotations/annotations_genes_QuickGO.txt", 'r') as f1:
    for i in f1.readlines():
        i = i.strip('\n')
        i = i.split('\t')
        MGI_ID = i[0]
        annot = i[1]
        if MGI_ID in annotation.keys():
            if annot not in annotation[MGI_ID]:
                annotation[MGI_ID].append(annot)

        else:
            annotation[MGI_ID] = [annot]

output = open("../../data/mouse/annotations/annotations_genes_QuickGO_uniq.txt", 'w')
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

"""
new_output = open("../../data/human/annotations/Gene_developmental_process_annotated.txt", 'w')
with open("../../data/human/annotations/Gene_developmental_process.txt", 'r') as f1:
    for line in f1.readlines():
        line = line.strip('\n')
        i = line.split('\t')
        MGI_ID = i[0]
        new_output.write(line + '\t' + str(annotation[MGI_ID]) + '\n')

new_output.close()
"""
