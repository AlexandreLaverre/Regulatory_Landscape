#!/bin/bash

export sp=$1
export sample=$2
export cluster=$3

##########################################################################

if [ ${cluster} = "pbil" ]; then
    export pathHiCup=/beegfs/data/${USER}/Regulatory_landscape/result/HiCup
    export pathFinalData=/beegfs/data/necsulea/RegulatoryLandscapesManuscript
    export pathFragmentCoords=${pathFinalData}/SupplementaryDataset1/${sp}
fi

if [ ${sp} = "human" ]; then
    export genome="hg38"
fi

if [ ${sp} = "mouse" ]; then
    export genome="mm10"
fi

##########################################################################

if [ -e ${pathHiCup}/frag_coords_${genome}_formatted.bed ]; then
    echo "fragment coords already done"
else
    ## file format: remove ^chr from chromsome names
    cat ${pathFragmentCoords}/frag_coords_${genome}.bed | sed '1d' |  sed -e 's/^chr//'> ${pathHiCup}/frag_coords_${genome}_formatted.bed
fi

##########################################################################

for file in `ls ${pathHiCup}/${sample} | grep sorted.dedup.bam`
do
    export prefix=`basename ${file} .sorted.dedup.bam`
    bedtools coverage -a ${pathHiCup}/frag_coords_${genome}_formatted.bed -b ${pathHiCup}/${sample}/${file} -split > ${pathHiCup}/${sample}/${prefix}_frag_coverage.txt 
done

##########################################################################

