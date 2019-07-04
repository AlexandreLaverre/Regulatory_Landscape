#!/usr/bin/env python3
# coding=utf-8


import random

ref = "mouse"
sp = "opossum"
path = "/home/laverre/Documents/Regulatory_Landscape/"

output = open(path+"result/conservation/interfrag_distance_" + ref + "2" + sp + "_PECAN_echant.txt", 'w')
with open(path+"result/conservation/interfrag_distance_" + ref + "2" + sp + "_PECAN.txt", 'r') as f1:
    line = f1.readline()
    while line:
        nb = random.randint(1, 23)
        if nb > 10:
            output.write(line + '\n')

        line = f1.readline().strip('\n')

output.close()

