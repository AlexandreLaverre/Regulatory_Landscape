#!/usr/bin/env python3
# coding=utf-8


import random

ref = "mouse"
sp = "human"
path = "/home/laverre/Documents/Regulatory_Landscape/"

output = open(path+"result/conservation/interfrag_distance_" + ref + "2" + sp + "_echant2.txt", 'w')
with open(path+"result/conservation/interfrag_distance_" + ref + "2" + sp + ".txt", 'r') as f1:
    line = f1.readline()
    while line:
        nb = random.randint(1, 510000000)
        if nb < 1000000:
            output.write(line + '\n')

        line = f1.readline().strip('\n')

output.close()

