#!/usr/bin/env python3
# coding=utf-8

from scipy.stats import wasserstein_distance

f2 = open("../data/mouse/dist_real.txt", 'r')
real = []
for i in f2.readlines():
    real.append(i)

f1 = open("../data/mouse/dist_simul.txt", 'r')
simul = []
for i in f1.readlines():
    simul.append(i)

print(wasserstein_distance(real, simul))



