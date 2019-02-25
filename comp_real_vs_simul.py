#!/usr/bin/env python3
# coding=utf-8

from matplotlib import pyplot as plt
import numpy as np

sp = "mouse" ##sys.argv[1]

contact = {}
frag_size = []
vect_dist = []
trans = 0
too_large = 0
too_short = 0
total = 0
with open("../data/"+sp+"/all_interactions/all_interactions.txt") as f3:
    for i in f3.readlines()[1:]:
        i = i.split("\t")
        midbait = ((int(i[1]) + int(i[2])) / 2)
        midcontact = ((int(i[4]) + int(i[5])) / 2)
        frag = (str(i[0]) + ":" + str(i[1]) + ":" + str(i[2]))
        PIR = (str(i[4]) + ":" + str(i[5]))
        total += 1

        if i[0] == i[3]:
            if 25000 < abs(midbait-midcontact) < 10000000:
                vect_dist.append(midbait - midcontact)

                if frag in contact.keys():
                    contact[frag].append(PIR)
                else:
                    contact[frag] = [PIR]
                    frag_size.append(int(i[2]) - int(i[1]))
            if 25000 > abs(midbait-midcontact):
                too_short += 1
            if abs(midbait-midcontact) > 10000000:
                too_large += 1
        else:
            trans += 1

nb_contact = []
for i in contact.keys():
    nb_contact.append(len(contact[i]))

print("###### Real ######")
print("Nombre d'interactions total:", total)
print("Nombre d'interactions valides", len(vect_dist))
print("Too large:", too_large, "soit", (too_large/total)*100, "%")
print("Too short:", too_short, "soit", (too_short/total)*100, "%")
print("Trans:", trans, "soit", (trans/total)*100, "%")
print()
print("Distance moyenne", np.mean(vect_dist))
print("Distance mediane", np.median(vect_dist))
plt.hist(vect_dist, min(len(set(vect_dist)), 500), (-2000000, 2000000))
plt.title('Distance distribution')
plt.xlabel('nucleotidic distance')
plt.ylabel('Frequency')
plt.show()

print("Frag size moyenne", np.mean(frag_size))
print("Frag size mediane", np.median(frag_size))
print("Frag size max", np.max(frag_size))
plt.hist(frag_size, min(len(set(frag_size)), 500), (0, 40000))
plt.title('Frag size distribution')
plt.xlabel('Frag size (pb)')
plt.ylabel('Frequency')
plt.show()

print("Nombre de contact moyen", np.mean(nb_contact))
print("Nombre de contact median", np.median(nb_contact))
plt.hist(nb_contact, min(len(set(nb_contact)), 200), (0, 200))
plt.title('Contact nb distribution')
plt.xlabel('Contact number')
plt.ylabel('Frequency')
plt.show()

contact_simul = {}
frag_size_simul = []
vect_dist_simul = []
midbait_out = 0
with open("../data/"+sp+"/Simulations/simulations_mouse_10Mb_bin5kb_fragoverbin.txt") as f3:
    for i in f3.readlines()[1:]:
        i = i.strip("\n")
        i = i.split("\t")
        midbait = ((int(i[1]) + int(i[2])) / 2)
        midcontact = ((int(i[3]) + int(i[4])) / 2)
        frag = (str(i[0]) + ":" + str(i[1]) + ":" + str(i[2]))
        PIR = (str(i[3]) + ":" + str(i[4]))

        if 25000 < abs(midbait - midcontact) < 10000000:
            vect_dist_simul.append(midbait - midcontact)

            if frag in contact_simul.keys():
                contact_simul[frag].append(PIR)
            else:
                contact_simul[frag] = [PIR]
                frag_size_simul.append(int(i[2]) - int(i[1]))
        else:
            midbait_out += 1

nb_contact_simul = []
for i in contact_simul.keys():
    nb_contact_simul.append(len(contact_simul[i]))

print("###### Simulations ######")
print("Nombre d'interactions", len(vect_dist_simul))
print("Distance moyenne", np.mean(vect_dist_simul))
print("Distance mediane", np.median(vect_dist_simul))
plt.hist(vect_dist_simul, min(len(set(vect_dist_simul)), 500), (-2000000, 2000000))
plt.title('Distance distribution')
plt.xlabel('nucleotidic distance')
plt.ylabel('Frequency')
plt.show()

print("Frag size moyenne", np.mean(frag_size_simul))
print("Frag size mediane", np.median(frag_size_simul))
print("Frag size max", np.max(frag_size_simul))
plt.hist(frag_size_simul, min(len(set(frag_size_simul)), 500), (0, 40000))
plt.title('Frag size distribution')
plt.xlabel('Frag size (pb)')
plt.ylabel('Frequency')
plt.show()

print("Nombre de contact moyen", np.mean(nb_contact_simul))
print("Nombre de contact median", np.median(nb_contact_simul))
plt.hist(nb_contact_simul, min(len(set(nb_contact_simul)), 200), (0, 200))
plt.title('Contact nb distribution')
plt.xlabel('Contact number')
plt.ylabel('Frequency')
plt.show()

print("Nb out of range:", midbait_out)


# Similitude
frag_absent = []
same_interaction = []
for frag in contact.keys():
    if frag not in contact_simul.keys():
        frag_absent.append(frag)
    else:
        for PIR in contact[frag]:
            if PIR in contact_simul[frag]:
                same_interaction.append(str(frag)+':'+str(PIR))


print("###### Similitudes #####")
print("Nb frag real", len(contact.keys()))
print("Nb frag simul", len(contact_simul.keys()))
print("Nb frag absent de la simul:", len(frag_absent))
print("Frag absent de la simul:", frag_absent)
print()
print("Nb total interactions:", sum(len(contact[i]) for i in contact.keys()))
print("Nb total interactions simulées:", sum(len(contact_simul[i]) for i in contact_simul.keys()))
print("Nb total interactions simulées uniques:", sum(len(set(contact_simul[i])) for i in contact_simul.keys()), "soit", (sum(len(set(contact_simul[i])) for i in contact_simul.keys())/sum(len(contact_simul[i]) for i in contact_simul.keys()))*100, "%")
print("Nb interactions absentes:", sum(len(contact[i]) for i in contact.keys())-sum(len(contact_simul[i]) for i in contact_simul.keys()), "soit", (((sum(len(contact[i]) for i in contact.keys())-sum(len(contact_simul[i]) for i in contact_simul.keys()))/sum(len(contact_simul[i]) for i in contact_simul.keys())))*100, "%")
print("Nb interactions identiques:", len(same_interaction), "soit", (len(same_interaction)/sum(len(contact_simul[i]) for i in contact_simul.keys()))*100, "%")

#print("Nb total interactions simulées:", sum(len(contact_simul[i]) for i in contact_simul.keys()))

dist_real = open("../data/"+sp+"/dist_real.txt", "w")
dist_simul = open("../data/"+sp+"/Simulations/dist_simulations_mouse_10Mb_bin5kb_fragoverbin.txt", "w")
# Output comparaison distribution dans R
for i in vect_dist:
    dist_real.write(str(i)+'\n')

for i in vect_dist_simul:
    dist_simul.write(str(i)+'\n')

dist_real.close()
dist_simul.close()
