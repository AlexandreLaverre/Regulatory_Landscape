#!/usr/bin/env python
# coding=utf-8

from Bio import SeqIO
import os
import sys

file = sys.argv[1]
sp1 = sys.argv[2]
sp2 = sys.argv[3]

path = "/mnt/result/genome_alignment/"+sp1+"2other/"
sp1_genome = path+"mouse_restriction_fragments.fa"
sp1_dict = SeqIO.to_dict(SeqIO.parse(sp1_genome, "fasta"))

sp2_genome = path+sp1+"2"+sp2+"_restriction_fragments_0.1.fa"
sp2_dict = SeqIO.to_dict(SeqIO.parse(sp2_genome, "fasta"))

error_file = open(path+file+"_frag_error.txt", "w")

with open(path+"after/"+file, 'r') as f1:
    for i in f1.readlines():
        i = i.strip("\n")
        i = i.split("\t")

        os.mkdir(path+'sp1_frag_%s/' % i[1])
        os.chdir(path+'sp1_frag_%s/' % i[1])

        if str(i[1]) in sp1_dict.keys() and str(i[1]) in sp2_dict.keys():
            sp1_seq = sp1_dict[i[1]]
            sp1_seq.id = str(sp1)
            sp1_seq.description = ""
            output_sp1 = open("./"+sp1 , "w+")
            SeqIO.write(sp1_seq, output_sp1, "fasta")
            output_sp1.close()

            sp2_seq = sp2_dict[i[1]]
            sp2_seq.id = str(sp2)
            sp2_seq.description = ""
            output_sp2 = open("./"+sp2 , "w+")
            SeqIO.write(sp2_seq, output_sp2, "fasta")
            output_sp2.close()

            os.system("all_bz '("+sp1+" "+sp2+")' ")
            os.system("tba '("+sp1+" "+sp2+")' *.*.maf tba.maf")
            os.system("gzip tba.maf")
            os.system("mv tba.maf.gz "+path+"tba_alignments/"+sp1+"2"+sp2+"/"+i[1]+"."+i[0]+".maf.gz")
            os.system("rm -r "+path+"sp1_frag_"+i[1])

        else:
            error_file.write(str(i)+"\n")


error_file.close()


