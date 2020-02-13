#!/bin/bash

################################################################################################################################################

specie=$1
path=/home/laverre/Documents/Regulatory_Landscape
pathScripts=${path}/scripts/overlap_genomic_elements
pathOverlap=${path}/data/${specie}/overlap
bonus=$2

if [ ${specie} = "human" ]; then
enhancers_files=(potential_enhancers/RoadMap_enhancer_genomic_positions_hg38.bed potential_enhancers/GRO_seq_enhancer_genomic_positions_hg38.bed potential_enhancers/ENCODE_enhancer_genomic_positions_hg38.bed potential_enhancers/CAGE_enhancer_genomic_positions_hg38.bed)
converted_frag=mouse2human_merged_interacted_fragments_bait_corrected.bed
converted_prefix=mouse2human_merged

else
enhancers_files=(potential_enhancers/ENCODE_enhancers_genomic_positions_mm10.bed potential_enhancers/CAGE/CAGE_enhancers_genomic_positions_mm10.bed)
converted_frag=human2mouse_merged_interacted_fragments_bait_corrected.bed
converted_prefix=human2mouse_merged
fi

reference_file=${specie}_merged_interacted_fragments_bait_corrected.txt
prefix=${specie}_merged
interest_files=(annotations/coding_exons.bed annotations/nocoding_exons.bed annotations/all_exons.bed annotations/genes.bed annotations/TSS_genes.bed annotations/phastconselements_noexonic250.bed annotations/repeatmasker.bed)

################################################################################################################################################

# Overlap with interest files
for file in "${interest_files[@]}"
do
suffix=${file%.*}
suffix=${suffix##*/}
if test -f "${pathOverlap}/${prefix}_overlap_${suffix}.txt"; then
echo "############ Overlap between ${prefix} vs ${suffix} already done ! ############"
else
echo "############ Running ${prefix} vs ${suffix} ############"
${pathScripts}/overlap.py ${specie} ${reference_file} ${file} overlap/${prefix}_overlap_${suffix}.txt --intraoverlap --countbp -v
fi
done

# Overlap with enhancers
for enh in "${enhancers_files[@]}"
do
suffix=${enh%_enhancer*}
suffix=${suffix##*/}
if test -f "${pathOverlap}/${prefix}_overlap_${suffix}.txt"; then
echo "############ Overlap between ${prefix} vs ${suffix} already done ! ############"
else
echo "############ Running ${prefix} vs ${suffix} ############"
${pathScripts}/overlap.py ${specie} ${reference_file} ${enh} overlap/${prefix}_overlap_${suffix}.txt --intraoverlap --countbp -v
fi
done

# Overlap lifted fragments to reference species
suffix=${reference_file%_interacted*}
if test -f "${pathOverlap}/${converted_prefix}_overlap_${suffix}.txt"; then
echo "############ Overlap between ${converted_prefix} vs ${suffix} already done ! ############"
else
echo "############ Running ${converted_prefix} vs ${suffix} ############"
${pathScripts}/overlap.py ${specie} ${converted_frag} ${reference_file} overlap/${converted_prefix}_overlap_${suffix}.txt -v
fi

if [ ${bonus} = "T" ]; then
${pathScripts}/overlap.py ${specie} ${reference_file} annotations/all_exons.bed overlap/${prefix}_overlap_all_exons250.txt --extend 250 --intraoverlap --countbp -v
${pathScripts}/overlap.py ${specie} ${reference_file} annotations/TSS_genes.bed overlap/${prefix}_overlap_TSS_1Kb.txt --extend 1000 --intraoverlap --countbp -v
${pathScripts}/overlap.py ${specie} ${reference_file} annotations/genes.bed overlap/${prefix}_overlap_genes_1Kb.txt --extend 1000 --intraoverlap --countbp -v
fi



