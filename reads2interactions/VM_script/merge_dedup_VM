#!/bin/bash
specie=$1			#ex : mouse
author_ctype=$2			#ex : Schoenfelder/ESC
name=$3				#ex : ESC_rep1
path_result="mnt/data/${specie}/${author_ctype}/"

samtools merge -n ${path_result}${name}_sorted.bam ${path_result}SRR5799157/SRR5799157_1_2.hicup.bam ${path_result}SRR5799158/SRR5799158_1_2.hicup.bam ${path_result}SRR5799159/SRR5799159_1_2.hicup.bam -@ 10

hicup_deduplicator --zip ${path_result}${name}_sorted.bam  --outdir ${path_result} --threads 10



