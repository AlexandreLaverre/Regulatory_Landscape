#!/bin/bash
specie=$1
genome=$2
author_ctype=$3

path_bam="/mnt/data/${specie}/${author_ctype}/"
path_data="/mnt/data/CHICAGO_files/${specie}/${genome}/"
path_result="/mnt/result/CHICAGO/${specie}/${author_ctype}/"
mkdir -p ${path_result}
for sample in `ls ${path_bam}*.dedup.bam`
do
echant=$(echo ${sample} | rev | cut -d / -f 1  | rev | cut -d . -f 1)
bam2chicago.sh ${sample} ${path_data}Digest_${genome}_HindIII_None.txt.baitmap ${path_data}Digest_${genome}_HindIII_None.txt.rmap ${path_result}${echant}
done

