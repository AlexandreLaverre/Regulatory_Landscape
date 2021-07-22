#!/bin/bash

################################################################################

species=$1
way=$2

enhancers=(FANTOM5 ENCODE restriction_fragments)

if [ "${species}" = "human" ]; then
enhancers=("${enhancers[@]}" RoadmapEpigenomics FOCS_GRO_seq)
fi

################################################################################


for seq in ${enhancers[@]}
do
	path=/beegfs/data/alaverre/Regulatory_landscape/result/phyloP/${species}/${seq}
	echo "Formatting ${species} ${seq}..."

	if [ -e ${path}/phyloP_${way}_MaskedExons_Ensembl94.txt ]; then
		echo "already done"
	else
		grep -hv CoveredLength ${path}/*_${way}_chr* > ${path}/phyloP_${way}_MaskedExons_Ensembl94.txt
		cut -f 2,3,4 ${path}/phyloP_${way}_MaskedExons_Ensembl94.txt > ${path}/ID
		cut -f 2-7 ${path}/phyloP_${way}_MaskedExons_Ensembl94.txt > ${path}/other
		sed -i "s/\t/:/g" ${path}/ID
		paste -d "\t" ${path}/ID ${path}/other > ${path}/phyloP_${way}_MaskedExons_Ensembl94.txt
		sed -i '1 i\ID\tChr\tStart\tEnd\tScore\tCoveredLength\tAnalyzedLength' ${path}/phyloP_${way}_MaskedExons_Ensembl94.txt
		rm ${path}/ID ${path}/other ${path}/*_${way}_chr*
	fi
done

################################################################################
