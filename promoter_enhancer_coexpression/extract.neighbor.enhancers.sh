#!/bin/bash

export sp=$1
export dataset=$2
export cluster=$3

###################################################################################

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/necsulea/RegulatoryLandscapes
    export pathManuscript=/beegfs/data/necsulea/RegulatoryLandscapesManuscript
fi

export pathFOCS=${path}/data/FOCS/${sp}/${dataset}
export pathOriginalInteractions=${pathManuscript}/SupplementaryDataset1/${sp}
export pathSimulatedInteractions=${pathManuscript}/SupplementaryDataset2/${sp}
export pathResults=${path}/results/co_expression_analysis/${sp}/${dataset}
export pathScripts=${path}/scripts/promoter_enhancer_coexpression

###################################################################################

if [ -e ${pathResults} ]; then
    echo "output path already there"
else
    mkdir -p ${pathResults} 
fi

###################################################################################

if [ ${sp} = "human" ]; then
    export assembly="hg38"
fi

if [ ${sp} = "mouse" ]; then
    export assembly="mm10"
fi

###################################################################################

if [ -e ${pathResults}/neighbor_promoters_enhancers.txt ]; then
    echo "already done"
else
    perl ${pathScripts}/extract.neighbor.enhancers.pl  --pathPromoters=${pathFOCS}/promoter_coordinates_${assembly}.bed  --pathEnhancers=${pathFOCS}/enhancer_coordinates_${assembly}.bed --maxDistance=1000000 --pathOutput=${pathResults}/neighbor_promoters_enhancers.txt 
fi

###################################################################################
