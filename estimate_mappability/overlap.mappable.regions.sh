#!/bin/bash

export sp=$1
export readlen=$2
export cluster=$3

#################################################################################

if [ ${cluster} = "pbil" ]; then
    export pathFinalData=/beegfs/data/necsulea/RegulatoryLandscapesManuscript
    export path=/beegfs/data/${USER}/RegulatoryLandscapes
fi

export pathFragments=${pathFinalData}/SupplementaryDataset1/${sp}
export pathResults=${path}/results/mappability/readlength${readlen}/${sp}
export pathScripts=${path}/scripts/estimate_mappability

#################################################################################

if [ ${sp} = "human" ]; then
    export genome="hg38"
fi

if [ ${sp} = "mouse" ]; then
    export genome="mm10"
fi

#################################################################################

perl ${pathScripts}/overlap.mappable.regions.pl --pathRestrictionFragments=${pathFragments}/frag_coords_${genome}.bed --pathMappedRegions=${pathResults}/mapped_regions.txt --margin=500 --step=10 --pathOutput=${pathResults}/overlap_restriction_fragments_mapped_regions_margin500_step10.txt

#################################################################################
