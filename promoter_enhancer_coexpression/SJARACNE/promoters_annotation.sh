#!/bin/bash

##############################################################################################

species=$1 # human or mouse
enh=$2 # FANTOM5 RoadmapEpigenomics ENCODE GRO-seq

path=/beegfs/data/alaverre/Regulatory_landscape/

pathData=/beegfs/data/necsulea/RegulatoryLandscapes/data/FOCS/${species}/${enh}
pathAnnot=${path}/result/Supplementary_dataset3_annotations/${species}/human_TSS_Ensembl94.bed
pathOutput=${path}/data/activities_correlation/${enh}/

pathScripts=${path}/scripts/overlap_genomic_elements

############################################################################################

echo "Running overlap for : ${enh} in ${species}"

${pathScripts}/overlap.py ${pathData}/promoter_coordinates_hg38.bed ${pathAnnot} ${pathOutput}/promoter_hg38_annotated_TSS_1kb.bed --reference_ID --interest_ID --extend 1000 -v

############################################################################################

echo "Get filtered promoters activity"

grep -v "NA" ${pathOutput}/promoter_hg38_annotated_TSS_1kb.bed | cut -f 1 > ${pathOutput}/promoter_annotated_TSS_1kb_filtered.bed
grep -Fwf ${pathOutput}/promoter_annotated_TSS_1kb_filtered.bed ${pathData}/promoter_activity.txt > ${pathOutput}/promoter_activity_filtered.txt


############################################################################################

echo "Formatting activities files for SJARACNE"

# Concatenate enhancers and promoters activity
cat ${pathData}/enhancer_activity.txt ${pathOutput}/promoter_activity_filtered.txt > ${pathOutput}/enh_prom_activity_SJARACNE.txt
sed -i '1d' ${pathOutput}/enh_prom_activity_SJARACNE.txt

# Add isoformID column
tail -n +2 ${pathData}/enhancer_activity.txt | cut -f 1  > enh_ID
cut -f 1 ${pathOutput}/promoter_activity_filtered.txt > ${pathOutput}/prom_ID.txt

if [ ${enh} == "FANTOM5" ]; then
sed -i "s/-/_/g" enh_ID 
sed -i "s/\../_/g" ${pathOutput}/prom_ID.txt

elif [ ${enh} == "ENCODE" ] || [ ${enh} == "RoadmapEpigenomics" ] ; then
sed -i "s/e_/e./g" enh_ID 
sed -i "s/p_/p./g" ${pathOutput}/prom_ID.txt 

elif [ ${enh} == "GRO-seq" ]; then
sed -i "s/e_/e./g" enh_ID 
sed -i "s/\n/_1\n/g" ${pathOutput}/prom_ID.txt 
fi

cat enh_ID ${pathOutput}/prom_ID.txt > all_ID
paste all_ID ${pathOutput}/enh_prom_activity_SJARACNE.txt > ${pathOutput}/enh_prom_activity_SJARACNE.txt2

# Add sample row
head -n 1 ${pathData}/promoter_activity.txt > samples
sed -i "s/^/isoformId\tgeneSymbol\t/" samples
cat samples ${pathOutput}/enh_prom_activity_SJARACNE.txt2 > ${pathOutput}/enh_prom_activity_SJARACNE.txt

rm samples ${pathOutput}/enh_prom_activity_SJARACNE.txt2




