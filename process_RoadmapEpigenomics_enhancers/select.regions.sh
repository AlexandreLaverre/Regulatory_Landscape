#!/bin/bash

export type=$1
export cluster=$2

###################################################################################

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/necsulea/RegulatoryLandscapes
else
    if [ ${cluster} = "in2p3" ]; then
	export path=/sps/biometr/necsulea/RegulatoryLandscapes
    else
	if [ ${cluster} = "in2p3_local" ]; then
	    export path=/sps/biometr/necsulea/RegulatoryLandscapes
	else
	    echo "unknown cluster"
	    exit
	fi
    fi
fi

export pathRoadmap=${path}/data/RoadmapEpigenomics/hg19
export pathCoverage=${pathRoadmap}/signal_tracks
export pathEnhancers=${pathRoadmap}/regulatory_regions/enhancer_regions
export pathPromoters=${pathRoadmap}/regulatory_regions/promoter_regions
export pathScripts=${path}/scripts/process_RoadmapEpigenomics_enhancers

if [ ${type} = "enhancers" ]; then
    export pathResults=${pathEnhancers}
    export prefix="all_enhancers"
else
    if [ ${type} = "promoters" ]; then
	export pathResults=${pathPromoters}
	export prefix="all_promoters"
    else
	echo "unknown type"
	exit
    fi
fi

###################################################################################

perl ${pathScripts}/select.regions.pl --pathCombinedCoverage=${pathResults}/${prefix}_coverage_allsamples.txt --minCoverage=5 --minSamples=5 --pathOutput=${pathResults}/selected_${type}_minFC5_minSamples5_coverage_allsamples.bed

###################################################################################
