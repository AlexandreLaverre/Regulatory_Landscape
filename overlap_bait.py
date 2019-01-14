#!/usr/bin/env python3
# coding=utf-8

n = 1
not_chevauch = 0
chevauch = []

with open("../data/Digest_mm9_HindIII_baits_ID.baitmap", 'r') as f1:
    for i in f1.readlines()[1:]:
        i = i.strip('\n')
        i = i.split("\t")

        if n == 1:
            previous_chr = i[0]
            previous_start = int(i[1])
            previous_end = int(i[2])

        else:
            current_chr = i[0]
            current_start = int(i[1])
            current_end = int(i[2])

            if current_chr == previous_chr:
                if current_start >= previous_end:
                    not_chevauch += 1
                else:
                    chevauch.append(current_chr+'\t'+current_start+'\t'+current_end)
            else:
                not_chevauch += 1

            current_chr = previous_chr
            current_start = previous_start
            current_end = previous_end
        n += 1

print("Nb de bait:", n)
print("Nb de bait overlap:", len(chevauch))
print("Nb de bait ok:", not_chevauch)
