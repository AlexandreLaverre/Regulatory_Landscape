#!/bin/bash

################################################################################################################################################

specie=$1
path=/home/laverre/Documents/Regulatory_Landscape
pathScripts=${path}/scripts/overlap_genomic_elements
pathOverlap=${path}/data/${specie}/overlap
data=$2
bonus=$3

if [ ${specie} = "human" ]; then
enhancers_files=(potential_enhancers/RoadMap_enhancer_genomic_positions_hg38.bed potential_enhancers/GRO_seq_enhancer_genomic_positions_hg38.bed potential_enhancers/ENCODE_enhancer_genomic_positions_hg38.bed potential_enhancers/CAGE_enhancer_genomic_positions_hg38.bed)
converted_enhancers=(potential_enhancers/mouse2human_ENCODE.bed potential_enhancers/mouse2human_CAGE.bed)
if [ ${data} = "merged" ]; then
converted_frag=mouse2human_corrected_merged.bed
converted_prefix=mouse2human_corrected_merged
elif [ ${data} = "samples" ]; then
converted_frag=mouse2human_samples_simulations.bed
converted_prefix=mouse2human_samples_simulations
else
converted_frag=mouse2human_lifted_frag.bed
converted_prefix=mouse2human_lifted_frag
fi

else
enhancers_files=(potential_enhancers/ENCODE_enhancers_genomic_positions_mm10_extend.bed potential_enhancers/CAGE/CAGE_enhancers_genomic_positions_mm10.bed)
converted_enhancers=(potential_enhancers/human2mouse_ENCODE.bed potential_enhancers/human2mouse_CAGE.bed potential_enhancers/human2mouse_GRO_seq.bed potential_enhancers/human2mouse_RoadMap.bed)
if [ ${data} = "merged" ]; then
converted_frag=human2mouse_corrected_merged.bed
converted_prefix=human2mouse_corrected_merged
elif [ ${data} = "samples" ]; then
converted_frag=human2mouse_samples_simulations.bed
converted_prefix=human2mouse_samples_simulations
else
converted_frag=human2mouse_lifted_frag.bed
converted_prefix=human2mouse_lifted_frag
fi
fi

if [ ${data} = "merged" ]; then
reference_file=${specie}_corrected_merged.txt
prefix=${specie}_corrected_merged
elif [ ${data} = "samples" ]; then
reference_file=${specie}_samples_simulations.txt
prefix=${specie}_samples_simulations
else
reference_file=${specie}_frag_coord.txt
prefix=${specie}_all_fragments
fi

interest_files=(annotations/coding_exons.bed annotations/nocoding_exons.bed annotations/all_exons.bed annotations/genes.bed annotations/TSS_genes.bed annotations/phastconselements_noexonic250.bed annotations/repeatmasker.bed)

################################################################################################################################################

# Overlap interest files
for file in "${interest_files[@]}"
do
suffix=${file%.*}
suffix=${suffix##*/}
if test -f "${pathOverlap}/${prefix}_overlap_${suffix}.txt"; then
echo "############ Overlap between ${prefix} vs ${suffix} already done ! ############"
else
echo "############ Running ${prefix} vs ${suffix} ############"
${pathScripts}/overlap.py ${specie} ${reference_file} ${file} overlap/${prefix}_overlap_${suffix}.txt --intraoverlap --count_overlap -v
fi
done

# Overlap enhancers
for enh in "${enhancers_files[@]}"
do
suffix=${enh%_enhancer*}
suffix=${suffix##*/}
if test -f "${pathOverlap}/${prefix}_overlap_${suffix}.txt"; then
echo "############ Overlap between ${prefix} vs ${suffix} already done ! ############"
else
echo "############ Running ${prefix} vs ${suffix} ############"
${pathScripts}/overlap.py ${specie} ${reference_file} ${enh} overlap/${prefix}_overlap_${suffix}.txt --count_overlap -v
fi
done

# Overlap enhancers with overlap
for enh in "${enhancers_files[@]}"
do
suffix=${enh%_enhancer*}
suffix=${suffix##*/}
if test -f "${pathOverlap}/${prefix}_overlap_${suffix}_intraoverlap.txt"; then
echo "############ Overlap between ${prefix} vs ${suffix} already done ! ############"
else
echo "############ Running ${prefix} vs ${suffix} intraoverlap ! ############"
${pathScripts}/overlap.py ${specie} ${reference_file} ${enh} overlap/${prefix}_overlap_${suffix}_intraoverlap.txt --intraoverlap --count_overlap -v
fi
done


# Overlap with converted enhancers
for enh in "${converted_enhancers[@]}"
do
suffix=${enh%.bed*}
suffix=${suffix##*/}
if test -f "${pathOverlap}/${prefix}_overlap_${suffix}.txt"; then
echo "############ Overlap between ${prefix} vs ${suffix} already done ! ############"
else
echo "############ Running ${prefix} vs ${suffix} ############"
${pathScripts}/overlap.py ${specie} ${reference_file} ${enh} overlap/${prefix}_overlap_${suffix}.txt --count_overlap -v
fi
done

# Overlap lifted fragments to reference species
suffix=${reference_file%_interacted*}
if test -f "${pathOverlap}/${converted_prefix}_overlap_${suffix}.txt"; then
echo "############ Overlap between ${converted_prefix} vs ${suffix} already done ! ############"
else
echo "############ Running ${converted_prefix} vs ${suffix} ############"
${pathScripts}/overlap.py ${specie} ${converted_frag} ${reference_file} overlap/${converted_prefix}_overlap_${suffix}.txt --count_overlap -v
fi

if [ ${bonus} = "T" ]; then
${pathScripts}/overlap.py ${specie} ${reference_file} annotations/all_exons.bed overlap/${prefix}_overlap_all_exons250.txt --extend 250 --intraoverlap --count_window -v
${pathScripts}/overlap.py ${specie} ${reference_file} annotations/TSS_genes.bed overlap/${prefix}_overlap_TSS_1Kb.txt --extend 1000 --intraoverlap --count_window -v
${pathScripts}/overlap.py ${specie} ${reference_file} annotations/genes.bed overlap/${prefix}_overlap_genes_1Kb.txt --extend 1000 --intraoverlap --count_window -v
fi



