#!/bin/bash
specie=$1
genome=$2
author_ctype=$3
path_data="/mnt/data/reads_PC_HIC/${specie}/${author_ctype}/"
path_result="/mnt/result/HiCup/${specie}/${author_ctype}/"

for sample in `ls ${path_data}*1.fastq.gz | xargs -n 1 basename | cut -d _ -f 1`
do
mkdir -p ${path_result}${sample}
## Config files for HiCUP ##
echo "Outdir: ${path_result}${sample}/" > ${path_result}${sample}/config_files.txt
echo "Zip: 1" >> ${path_result}${sample}/config_files.txt
echo "Keep: 0" >> ${path_result}${sample}/config_files.txt
echo "Threads: 10" >> ${path_result}${sample}/config_files.txt
echo "Bowtie2: /mnt/Tools/bowtie2/bowtie2" >> ${path_result}${sample}/config_files.txt 
echo "Digest: /mnt/data/genome_ref/Digest_${genome}_HindIII_None.txt" >> ${path_result}${sample}/config_files.txt
echo "Index: /mnt/data/genome_ref/${genome}" >> ${path_result}${sample}/config_files.txt
echo "R: /usr/bin/R" >> ${path_result}${sample}/config_files.txt
echo "Shortest: 150" >> ${path_result}${sample}/config_files.txt
echo "Longest: 800" >> ${path_result}${sample}/config_files.txt

echo "${path_data}${sample}_1.fastq.gz" >> ${path_result}${sample}/config_files.txt
echo "${path_data}${sample}_2.fastq.gz" >> ${path_result}${sample}/config_files.txt
echo "" >> ${path_result}${sample}/config_files.txt
sed -i '$d' ${path_result}${sample}/config_files.txt

## HiCUP ##
hicup --config  ${path_result}${sample}/config_files.txt

done
