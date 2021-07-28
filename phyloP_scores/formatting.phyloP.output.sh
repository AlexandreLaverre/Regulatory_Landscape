#!/bin/bash

################################################################################

species=$1
way=$2
masked=$3

enhancers=(FANTOM5 ENCODE restriction_fragments)

if [ "${species}" = "human" ]; then
enhancers=("${enhancers[@]}" RoadmapEpigenomics FOCS_GRO_seq)
fi

if [ "${masked}" = "TRUE" ]; then
export suffixExons=MaskedExons_Ensembl94
else
export suffixExons=Unmasked
fi

################################################################################


for seq in ${enhancers[@]}
do
	path=/beegfs/data/alaverre/Regulatory_landscape/result/phyloP/${species}/${seq}
	echo "Formatting ${species} ${seq}..."

	if [ -e ${path}/phyloP_${way}_${suffixExons}.txt ]; then
		echo "already done"
	else
		grep -hv CoveredLength ${path}/*_${way}_chr* > ${path}/phyloP_${way}_${suffixExons}.txt
		cut -f 2,3,4 ${path}/phyloP_${way}_${suffixExons}.txt > ${path}/ID
		cut -f 2-7 ${path}/phyloP_${way}_${suffixExons}.txt > ${path}/other
		sed -i "s/\t/:/g" ${path}/ID
		paste -d "\t" ${path}/ID ${path}/other > ${path}/phyloP_${way}_${suffixExons}.txt
		sed -i '1 i\ID\tChr\tStart\tEnd\tScore\tCoveredLength\tAnalyzedLength' ${path}/phyloP_${way}_${suffixExons}.txt
		rm ${path}/ID ${path}/other ${path}/*_${way}_chr*
	fi
done

################################################################################
