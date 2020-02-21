#! /usr/bin/env python

##### PIPELINE TO FIND AND ALIGN HOMOLOGOUS SEQUENCES BETWEEN SPECIES #####
# Data requirements:
# input file: data/sp_origin2other/seqs_coord.txt  list of interested coordinates
# soft masked genomes of interested species at data/genome_sm/sp.dna_sm.toplevel.fa.gz
# reference aligned genomes from USCS at data/sp2other/ref/sp1TOsp2.over.chain.gz
# exons coordinates of interested species at data/sp_all_exons.bed

# Scripts requirements:
# extract_frag_seq.py -> obtain DNA sequence from genomic coordinates
# overlap.py -> obtain the overlapping coordinates between 2 files (i.e: sequences and exons)
# pecan.py -> align pairs of homologous sequences
# extract.aln.stats.excluding.regions.pl -> obtain final alignment statistics

# Software requirements: LiftOver + PECAN + Snakemake

# Execute command : snakemake -j 200 --rerun-incomplete --cluster "sbatch -p normal -N 1 -o /beegfs/data/alaverre/Regulatory_landscape/test/result/mouse2other/mouse2human/log_files/slurm.out -e /beegfs/data/alaverre/Regulatory_landscape/test/result/mouse2other/mouse2human/log_files/slurm.err -c {params.threads} --mem={params.mem} -t {params.time}"

from snakemake.utils import listfiles

workdir: "/beegfs/data/alaverre/Regulatory_landscape/test/"
sp = "mouse"
interest_file = "merged_samples_simulated_fragments"
treshold = "0.8"


path="/beegfs/data/alaverre/Regulatory_landscape/test/"
pathScript=path+"script/"
result_path=path+"result/${sp}2other/duplication_rate/"
python_path = "/beegfs/data/soft/anaconda3/bin/python"

localrules : all, extract_seq, BLAT, extract_stats, merging

rule all:
    input :
        all_done

rule extract_seq:
    message:
        "Extracting masked sequences"
    output:
        "{result_path}/{sp}_restriction_fragments_mask/{sp}_{interest_file}.fa{sequence_file}"
    params:
        match = "0.1",
    shell:
        """
        {python_path} {pathScript}/extract_rm_seq_dupli.py {sp} {interest_file}
        """

rule BLAT:
    message:
        "BLAT : looking for duplicated sequences"
    input:
        "{result_path}/{sp}_restriction_fragments_mask/{sp}_{interest_file}.fa{sequence_file}",
        extract = result_path+"/log_files/extract_seq_"+sp+"_"+interest_file+"_done.txt"
    output:
        "{result_path}/output.psl/{sp}_{interest_file}.fa{sequence_file}_output.psl"
    params:
        time = "40:00:00", mem = "4G", threads = "1" ,
    shell:
        """
        {python_path} {pathScript}/run_blat.sh {sp} {interest_file}
        """

rule extract_stats:
    message:
        "Extract summary statistics of BLAT results"
    input:
       file = "{result_path}/output.psl/{sp}_{interest_file}.fa{sequence_file}_output.psl"
    output:
        "{result_path}/{sp}_{interest_file}.fa{sequence_file}_stats_{treshold}.txt"
    shell:
        """
        {python_path} {pathScript}/extract_duplication_stats.py {sp} {treshold} {input.file}
        """

def remap(wildcards):
    return listfiles(wildcards.sequence_file)

rule merging:
    message:
        "Make list of pairs ID to align"
    input:
        files =
        lifted = data_path+sp_origin+"2{sp_target}/"+sp_origin+"2{sp_target}_"+interest_file+".bed",
    output:
        ID_file = result_path+sp_origin+"2{sp_target}/"+sp_origin+"2{sp_target}_"+interest_file+"_ID.bed",
        done = result_path+sp_origin+"2{sp_target}/log_files/"+interest_file+"_ID_file_done"
    shell:
        """
        cat snakemake.utils.listfiles(
        """
