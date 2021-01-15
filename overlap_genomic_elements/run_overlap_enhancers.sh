#!/bin/bash

################################################################################################################################################

specie=$1
path=/home/laverre/Documents/Regulatory_Landscape
pathScripts=${path}/scripts/overlap_genomic_elements
pathOverlap=${path}/data/${specie}/overlap

if [ ${specie} = "human" ]; then
enhancers_files=(potential_enhancers/ENCODE_enhancer_genomic_positions_hg38.bed) #(potential_enhancers/RoadMap_enhancer_genomic_positions_hg38.bed potential_enhancers/GRO_seq_enhancer_genomic_positions_hg38.bed potential_enhancers/ENCODE_enhancer_genomic_positions_hg38.bed potential_enhancers/CAGE_enhancer_genomic_positions_hg38.bed)

else
enhancers_files=(potential_enhancers/ENCODE_enhancers_genomic_positions_mm10_extend.bed potential_enhancers/CAGE/CAGE_enhancers_genomic_positions_mm10.bed)
fi

interest_files=(annotations/all_exons.bed potential_enhancers/ENCODE_enhancer_genomic_positions_hg38.bed) #annotations/all_genes.bed

################################################################################################################################################

# Overlap interest files
#for file in "${interest_files[@]}"
#do
#suffix=${file%.*}
#suffix=${suffix##*/}
#	for enh in "${enhancers_files[@]}"
#	do
#	prefix=${enh%_enhancer*}
#	prefix=${prefix##*/}
#	if test -f "${pathOverlap}/${prefix}_overlap_${suffix}.txt"; then
#	echo "############ Overlap between ${prefix} vs ${suffix} already done ! ############"
#	else
#	echo "############ Running ${prefix} vs ${suffix} ############"
#	${pathScripts}/overlap.py ${specie} ${enh} ${file} overlap/${prefix}_overlap_${suffix}.txt --intraoverlap --count_overlap -v
#	fi
#	done
#done

# Overlap interest files +- distances
for file in "${interest_files[@]}"
do
suffix=${file%.*}
suffix=${suffix##*/}
	for enh in "${enhancers_files[@]}"
	do
	prefix=${enh%_enhancer*}
	prefix=${prefix##*/}
	if test -f "${pathOverlap}/${prefix}_overlap_${suffix}_50kb.txt"; then
	echo "############ Overlap between ${prefix} vs ${suffix} +- distances already done ! ############"
	else
	echo "############ Running ${prefix} vs ${suffix} +- distances ############"
	${pathScripts}/overlap.py ${specie} ${enh} ${file} overlap/${prefix}_overlap_${suffix}_50kb.txt --extend 50000 --intraoverlap --count_window -v
	${pathScripts}/overlap.py ${specie} ${enh} ${file} overlap/${prefix}_overlap_${suffix}_100kb.txt --extend 100000 --intraoverlap --count_window -v
	fi
	done
done




