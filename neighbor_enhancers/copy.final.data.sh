#!/bin/bash

############################################################

export path=/beegfs/data/necsulea/RegulatoryLandscapes
export pathFinalData=/beegfs/data/necsulea/RegulatoryLandscapesManuscript
export pathResults=${path}/results/neighbor_enhancers
export pathSuppDataset=${pathFinalData}/SupplementaryDataset4

export release=94

############################################################

for sp in human mouse
do
    cp ${pathResults}/${sp}/regulatory_regions_Ensembl${release}.txt ${pathSuppDataset}/${sp}/predicted_regulatory_regions_neighbor_TSS.txt
    
    for enh in `ls ${pathResults}/${sp}`
    do
	cp ${pathResults}/${sp}/${enh}/neighbor_enhancers_Ensembl${release}.txt ${pathSuppDataset}/${sp}/${enh}/predicted_enhancers_neighbor_TSS.txt
    done
done

############################################################


