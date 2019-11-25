#!/bin/bash

#####################################################################

export path=/sps/biometr/necsulea/RegulatoryLandscapes
export pathData=${path}/data/RoadmapEpigenomics
export pathAlignments=${path}/data/liftOver_files

#####################################################################

for type in promoter enhancer
do
    liftOver ${pathData}/hg19/${type}_regions/all_${type}s.bed ${pathAlignments}/hg19ToHg38.over.chain.gz ${pathData}/hg38/${type}_regions/all_${type}s.bed ${pathData}/hg38/${type}_regions/all_${type}s.unmapped
done

#####################################################################
