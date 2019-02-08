#!/usr/bin/env python3
# coding=utf-8

from matplotlib import pyplot as plt
import numpy as np

dic_dist = []
with open("/home/laverre/Documents/Regulatory_Landscape/data/mouse/all_interactions/all_interactions.txt", 'r') as f1:
    for i in f1.readlines()[1:]:
        i = i.split("\t")
        if i[0] == i[3]:
            midbait = ((int(i[1]) + int(i[2])) / 2)
            midcontact = ((int(i[4]) + int(i[5])) / 2)
            dist = midbait - midcontact
            if abs(dist) > 20000:
                dic_dist.append(dist)


print("Distance moyenne", np.mean(dic_dist))
print("Distance mediane", np.median(dic_dist))
print("Distance max", np.max(dic_dist))
print("Distance min", np.min(dic_dist))

plt.hist(dic_dist, min(len(set(dic_dist)), 500), (-1000000, 1000000))
plt.title('Distance distribution')
plt.xlabel('nucleotidic distance')
plt.ylabel('Frequency')
plt.show()
