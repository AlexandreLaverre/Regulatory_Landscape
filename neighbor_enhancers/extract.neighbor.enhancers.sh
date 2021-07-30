#!/bin/bash

export sp=$1
export dataset=$2
export cluster=$3

############################################################################

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/necsulea/RegulatoryLandscapes
    export pathFinalData=/beegfs/data/necsulea/RegulatoryLandscapesManuscript
fi

export pathScripts=${path}/scripts/neighbor_enhancers
export pathResults=${path}/results/neighbor_enhancers/${sp}
export pathEnhancers=${pathFinalData}/SupplementaryDataset4/${sp}/${dataset}

export release=94

############################################################################

perl ${pathScripts}/extract.neighbor.enhancers.pl --pathRegulatoryRegions=${pathResults}/regulatory_regions_Ensembl${release}.txt --pathEnhancers=${pathEnhancers}/enhancer_coordinates.bed --pathOutput=${pathResults}/${dataset}/neighbor_enhancers_Ensembl${release}.txt

############################################################################

