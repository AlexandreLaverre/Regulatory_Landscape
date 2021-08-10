#!/bin/bash

export sp=$1
export readlen=$2
export cluster=$3

##############################################################

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/${USER}/RegulatoryLandscapes
fi

export pathGenomeIndexes=${path}/data/genome_sequences/${sp}
export pathResults=${path}/results/mappability/readlength${readlen}/${sp}
export pathScripts=${path}/scripts/estimate_mappability

##############################################################

if [ -e ${pathResults}/mapped_regions.txt ]; then
    echo "result file already there, not doing anything"
    exit
fi

##############################################################

for chr in {1..32} X Y Z W 
do
    if [ -e ${pathResults}/fake_reads_chr${chr}_mappedregions.txt ]; then
	cat ${pathResults}/fake_reads_chr${chr}_mappedregions.txt >> ${pathResults}/mapped_regions.txt
    fi
done

##############################################################
