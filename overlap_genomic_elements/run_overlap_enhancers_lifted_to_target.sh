#!/bin/bash

################################################################################################################################################

specie=$1
path=/home/laverre/Documents/Regulatory_Landscape
pathScripts=${path}/scripts/overlap_genomic_elements
pathOverlap=${path}/data/${specie}/overlap

enhancers_files=(ENCODE CAGE)

if [ ${specie} = "human" ]; then
ref_sp="mouse"
genome="hg38"
else
ref_sp="human"
genome="mm10"
fi

################################################################################################################################################

# Overlap interest files

for enh in "${enhancers_files[@]}"
do
echo "############ Running ${enh} vs ${enh}_target ############"
${pathScripts}/overlap.py ${specie} potential_enhancers/${ref_sp}2${specie}_${enh}.bed potential_enhancers/${enh}_enhancer_genomic_positions_${genome}.bed overlap/${enh}_lifted_overlap_${enh}_target.txt -v --reference_ID
${pathScripts}/overlap.py ${specie} potential_enhancers/${ref_sp}2${specie}_${enh}.bed potential_enhancers/${enh}_enhancer_genomic_positions_${genome}.bed overlap/${enh}_lifted_overlap_${enh}_target_1kb.txt -v --extend 1000 --reference_ID
done




