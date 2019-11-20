#!/bin/bash

export cluster=$1

##############################################################

if [ ${cluster} = "in2p3" ] ; then
    export path=/sps/biometr/necsulea/RegulatoryLandscapes
fi
export pathRoadmap=${path}/data/RoadmapEpigenomics/hg19/regulatory_regions
export pathEnhancers=${pathRoadmap}/enhancer_regions
export pathPromoters=${pathRoadmap}/promoter_regions
export pathScripts=${path}/scripts/process_RoadmapEpigenomics_enhancers

##############################################################

export pathsE=""
export pathsP=""

for file in `ls ${pathEnhancers} | grep .bed.gz`
do
    export pathsE=${pathEnhancers}/${file},${pathsE}
done


for file in `ls ${pathPromoters} | grep .bed.gz`
do
    export pathsP=${pathPromoters}/${file},${pathsP}
done

echo ${pathsE}
echo ${pathsP}

##############################################################

perl ${pathScripts}/combine.coordinates.pl --pathsCoordinates=${pathsE} --pathOutput=${pathEnhancers}/all_enhancers.bed
perl ${pathScripts}/combine.coordinates.pl --pathsCoordinates=${pathsP} --pathOutput=${pathPromoters}/all_promoters.bed

##############################################################
