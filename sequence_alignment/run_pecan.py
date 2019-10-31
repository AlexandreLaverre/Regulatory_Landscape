#!/usr/bin/env python
# coding=utf-8

from Bio import SeqIO
import os
import sys

list = sys.argv[1] # ex: xaa
sp1 = sys.argv[2] # ex : human
sp2 = sys.argv[3] # ex : opossum
file = sys.argv[4] # ex : CAGE, merged_interacted_fragments,...

path = "/beegfs/data/alaverre/Regulatory_landscape/result/genome_alignment/"+sp1+"2other/"
sp1_genome = path+sp1+"_" + file +"_softmask.fa"
sp1_dict = SeqIO.to_dict(SeqIO.parse(sp1_genome, "fasta"))

path = "/beegfs/data/alaverre/Regulatory_landscape/result/genome_alignment/"+sp1+"2other/"+sp1+"2"+sp2+"/"
sp2_genome = path+sp1+"2"+sp2+ "_" + file +"_softmask.fa"
sp2_dict = SeqIO.to_dict(SeqIO.parse(sp2_genome, "fasta"))

error_file = open(path+'pecan_alignments/'+list+file+"_frag_error.txt", "w")
os.makedirs(path+'pecan_alignments/'+sp1+"2"+sp2 +"_" + file, exist_ok=True)

with open(path+"list/"+list, 'r') as f1: # Pairs seq aligned
    for i in f1.readlines()[1:]:
        i = i.strip("\n")
        i = i.split("\t")

        os.mkdir(path+'pecan_alignments/running/'+sp1+'_frag_%s/' % i[0])
        os.chdir(path+'pecan_alignments/running/'+sp1+'_frag_%s/' % i[0])

        if str(i[0]) in sp1_dict.keys() and str(i[0]) in sp2_dict.keys():
            sp1_seq = sp1_dict[i[0]]
            sp1_seq.id = str(sp1)
            sp1_seq.description = ""
            output_sp1 = open("./"+sp1, "w+")
            SeqIO.write(sp1_seq, output_sp1, "fasta")
            output_sp1.close()

            sp2_seq = sp2_dict[i[0]]
            sp2_seq.id = str(sp2)
            sp2_seq.description = ""
            output_sp2 = open("./"+sp2, "w+")
            SeqIO.write(sp2_seq, output_sp2, "fasta")
            output_sp2.close()

            os.system("java bp.pecan.Pecan -E '("+sp1+","+sp2+");' -F "+sp1+" "+sp2+" -G pecan.mfa")
            os.system("gzip pecan.mfa")
            os.system("mv pecan.mfa.gz " + path + "pecan_alignments/"+sp1+"2"+sp2 +"_" + file+"/" + i[0] + "." + i[1] + ".mfa.gz")
            os.system("rm -r "+path+'pecan_alignments/running/'+sp1+'_frag_'+i[0])

        else:
            error_file.write(str(i)+"\n")

error_file.close()
