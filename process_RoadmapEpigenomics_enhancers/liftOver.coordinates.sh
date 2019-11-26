#!/bin/bash

#####################################################################

export path=/sps/biometr/necsulea/RegulatoryLandscapes
export pathData=${path}/data/RoadmapEpigenomics
export pathAlignments=${path}/data/liftOver_files

#####################################################################

for type in promoter enhancer
do
    ## selected regions
    
    if [ -e  ${pathData}/hg38/regulatory_regions/${type}_regions/selected_${type}s.bed ]; then
	echo "already done"
    else
	liftOver ${pathData}/hg19/regulatory_regions/${type}_regions/selected_${type}s.bed ${pathAlignments}/hg19ToHg38.over.chain.gz ${pathData}/hg38/regulatory_regions/${type}_regions/selected_${type}s.bed ${pathData}/hg38/regulatory_regions/${type}_regions/selected_${type}s.unmapped
    fi

    ## all regions
    
    if [ -e  ${pathData}/hg38/regulatory_regions/${type}_regions/all_${type}s.bed ]; then
	echo "already done"
    else
	liftOver ${pathData}/hg19/regulatory_regions/${type}_regions/all_${type}s.bed ${pathAlignments}/hg19ToHg38.over.chain.gz ${pathData}/hg38/regulatory_regions/${type}_regions/all_${type}s.bed ${pathData}/hg38/regulatory_regions/${type}_regions/all_${type}s.unmapped
    fi

    
done

#####################################################################
