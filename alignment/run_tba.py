#!/usr/bin/env python3
# coding=utf-8

from Bio import SeqIO
import os
import sys

file=sys.argv[1]

path = "/beegfs/data/alaverre/result/genome_alignment/mouse2other/"
mouse_genome = path+"mouse_restriction_fragments.fa"
mouse_dict = SeqIO.to_dict(SeqIO.parse(mouse_genome, "fasta"))

human_genome = path+"mouse2human_restriction_fragments_0.1.fa"
human_dict = SeqIO.to_dict(SeqIO.parse(human_genome, "fasta"))

with open(path+file, 'r') as f1:
    for i in f1.readlines()[1:]:
        i = i.strip("\n")
        i = i.split("\t")
        os.mkdir(path+'mouse_frag_%s/' % i[0])
        os.chdir(path+'mouse_frag_%s/' % i[0])

        mouse_seq = mouse_dict[i[0]]
        mouse_seq.id = "mouse"
        mouse_seq.description = ""
        output_mouse = open("./mouse", "w+")
        SeqIO.write(mouse_seq, output_mouse, "fasta")
        output_mouse.close()

        human_seq = human_dict[i[0]]
        human_seq.id = "human"
        human_seq.description = ""
        output_human = open("./human" , "w+")
        SeqIO.write(human_seq, output_human, "fasta")
        output_human.close()

        os.system("all_bz '(mouse human)' ")
        os.system("tba '(mouse human)' *.*.maf tba.maf")
        os.system("gzip tba.maf")
        os.system("mv tba.maf.gz "+path+"tba_alignments/"+i[0]+"."+i[1]+".maf.gz")
        os.system("rm -r "+path+"mouse_frag_"+i[0])

