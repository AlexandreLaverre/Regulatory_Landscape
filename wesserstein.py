#!/usr/bin/env python3
# coding=utf-8

from scipy.stats import wasserstein_distance


def wasserstein(file):
    vect = []
    for i in file.readlines():
        vect.append(i)
    return vect


f1 = open("../data/mouse/Simulations/dist_real.txt", 'r')
real = wasserstein(f1)

f1 = open("../data/mouse/Simulations/dist_simul_simulations_5Mb_bin5kb_no_uniq.txt", 'r')
simul_no_uniq = wasserstein(f1)

f1 = open("../data/mouse/Simulations/dist_simul_simulations_5Mb_bin10kb_no_uniq.txt", 'r')
simul_no_uniq_10kb = wasserstein(f1)

f1 = open("../data/mouse/Simulations/dist_simul_simulations_2Mb_bin5kb_no_uniq.txt", 'r')
simul_no_uniq_2Mb = wasserstein(f1)

f1 = open("../data/mouse/Simulations/dist_simul_simulations_5Mb_bin5kb_uniq10_replace.txt", 'r')
simul_uniq10_replace = wasserstein(f1)

f1 = open("../data/mouse/Simulations/dist_simul_simulations_5Mb_bin5kb_uniq10_replace_prob.txt", 'r')
simul_uniq10_prob = wasserstein(f1)

f1 = open("../data/mouse/Simulations/dist_simul_alloverlap_uniq.txt", 'r')
simul_alloverlap = wasserstein(f1)

f1 = open("../data/mouse/Simulations/dist_simulations_mouse_10Mb_bin5kb_fragoverbin.txt", 'r')
simul_new = wasserstein(f1)

print("no_uniq : ", wasserstein_distance(real, simul_no_uniq))
print("no_uniq_10kb : ", wasserstein_distance(real, simul_no_uniq_10kb))
print("no_uniq_2Mb : ", wasserstein_distance(real, simul_no_uniq_2Mb))
print("uniq10_replace : ", wasserstein_distance(real, simul_uniq10_replace))
print("uniq10_prob : ", wasserstein_distance(real, simul_uniq10_prob))
print("alloverlap : ", wasserstein_distance(real, simul_alloverlap))
print("new : ", wasserstein_distance(real, simul_new))
f1.close()
