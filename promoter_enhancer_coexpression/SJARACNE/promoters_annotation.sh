#!/bin/bash

##############################################################################################

species=$1 # human or mouse
enh=$2 # FANTOM5 RoadmapEpigenomics ENCODE GRO-seq

path=/beegfs/data/alaverre/RegulatoryLandscapes/

pathData=/beegfs/data/necsulea/RegulatoryLandscapes/data/FOCS/${species}/${enh}
pathAnnot=${path}/result/Supplementary_dataset3_annotations/${species}/human_TSS_Ensembl94.bed
pathOutput=${path}/data/activities_correlation/${enh}/

pathScripts=${path}/scripts/overlap_genomic_elements

############################################################################################
echo "Running overlap for : ${enh} in ${species}"

${pathScripts}/overlap.py ${pathData}/promoter_coordinates_hg38.bed ${pathAnnot} ${pathOutput}/promoter_hg38_annotated_TSS_1kb.bed --reference_ID --interest_ID --extend 1000 -v

grep -v "NA" ${pathOutput}/promoter_hg38_annotated_TSS_1kb.bed | cut -f 1 > ${pathOutput}/promoter_annotated_TSS_1kb_filtered.bed
grep -Ff ${pathOutput}/promoter_annotated_TSS_1kb_filtered.bed ${pathData}/promoter_activity.txt > ${pathOutput}/promoter_activity_filtered.txt
