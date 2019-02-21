#!/usr/bin/env python3
# coding=utf-8

import os
from matplotlib import pyplot as plt

path = "/home/laverre/Documents/Regulatory_Landscape/data/mouse/CAGE/"

print("Calculating ENH length and tissu expression...")
size = {}
expr = {}
count_peak = {}
with open(path+"mouse.enhancers.expression.matrix", 'r') as f:
    first_line = f.readline()
    for i in f.readlines():
        i = i.strip("\n")
        i = i.split("\t")
        length = i[0].split(":")
        length = length[1].split("-")
        length = int(length[1]) - int(length[0])
        enh = i[0]
        size[enh] = length
        for tissu in range(1, len(i)):
            if tissu in expr.keys():
                expr[tissu] = float(expr[tissu]) + (int(i[tissu])/length)
            else:
                expr[tissu] = int(i[tissu])/length

            if int(i[tissu]) != 0:
                if tissu in count_peak.keys():
                    count_peak[tissu] = count_peak[tissu] + 1
                else:
                    count_peak[tissu] = 1


'''
count = 0
for i in count_peak.values():
    if i < 500:
        count += 1

print(count)

plt.hist(count_peak.values(), min(len(set(count_peak.values())), 50), (0, 2000))
plt.title('Tissu expression distribution')
plt.xlabel('TPM')
plt.ylabel('Frequency')
plt.show()
'''

print("Calculating TPM...")
tpm = {}
with open(path+"mouse.enhancers.expression.matrix", 'r') as f:
    for i in f.readlines()[1:]:
        i = i.strip("\n")
        i = i.split("\t")
        enh = i[0]
        tpm[enh] = []

        for tissu in range(1, len(i)):
            tpm[enh].append((int(i[tissu])/size[enh]) / expr[tissu])


print("Writting output...")
output = open(path+"mouse.enhancers.expression.tpm.matrix", 'w')
if os.stat(path+"mouse.enhancers.expression.tpm.matrix").st_size == 0:
    output.write(first_line)

for enh, exp in tpm.items():
    output.write(enh + "\t")
    count = 0
    for i in exp:
        count += 1
        if count == len(exp):
            output.write(str(i)+"\n")
        else:
            output.write(str(i) + "\t")

output.close()
