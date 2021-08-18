#!/bin/bash

################################################################################

species=$1
score=$2
way=$3
masked=$4

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
	path=/beegfs/data/alaverre/Regulatory_landscape/result/${score}/${species}/${seq}
	echo "Formatting ${species} ${seq}..."

	if [ -e ${path}/${score}_${way}_${suffixExons}.txt ]; then
		echo "already done"
	else
		grep -hv CoveredLength ${path}/*_${way}_chr* > ${path}/${score}_${way}_${suffixExons}.txt
		cut -f 2,3,4 ${path}/${score}_${way}_${suffixExons}.txt > ${path}/ID
		cut -f 2-7 ${path}/${score}_${way}_${suffixExons}.txt > ${path}/other
		sed -i "s/\t/:/g" ${path}/ID
		paste -d "\t" ${path}/ID ${path}/other > ${path}/${score}_${way}_${suffixExons}.txt
		sed -i '1 i\ID\tChr\tStart\tEnd\tScore\tCoveredLength\tAnalyzedLength' ${path}/${score}_${way}_${suffixExons}.txt
		rm ${path}/ID ${path}/other ${path}/*_${way}_chr*
	fi
done

################################################################################
