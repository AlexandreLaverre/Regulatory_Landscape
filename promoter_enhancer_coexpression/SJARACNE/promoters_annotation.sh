#!/bin/bash

##############################################################################################

species=$1 # human or mouse
enh=$2

path=/home/laverre/Documents/Regulatory_Landscape
pathScripts=${path}/scripts/overlap_genomic_elements
pathResults=FOCS_promoters/${enh}/

############################################################################################
echo "Running overlap for : ${enh} in ${species}"

${pathScripts}/overlap.py ${species} ${pathResults}promoter_coordinates_hg38.bed annotations/human_TSS_Ensembl94.bed ${pathResults}promoter_hg38_annotated_TSS_1kb.bed --reference_ID --interest_ID --extend 1000 -v
