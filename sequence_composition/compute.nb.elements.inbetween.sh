#!/bin/bash

export sp=$1

########################################################################

export pathOriginal=/beegfs/data/necsulea/RegulatoryLandscapesManuscript/SupplementaryDataset1/${sp}
export pathSimulated=/beegfs/data/necsulea/RegulatoryLandscapesManuscript/SupplementaryDataset2/${sp}
export pathAnnot=/beegfs/data/necsulea/RegulatoryLandscapesManuscript/SupplementaryDataset3/genes
export pathResults=/beegfs/data/necsulea/RegulatoryLandscapesManuscript/SupplementaryDataset5/${sp}
export pathScripts=/beegfs/data/necsulea/RegulatoryLandscapes/scripts/sequence_composition

########################################################################

export pathGenes=${pathAnnot}/${sp}_genes_Ensembl94.txt

########################################################################

for dataset in original simulated
do
    if [ ${dataset} = "original" ]; then
	export pathInteractions=${pathOriginal}/all_interactions.txt
    fi

    if [ ${dataset} = "simulated" ]; then
	export pathInteractions=${pathSimulated}/simulated_all_interactions.txt
    fi
    
    perl ${pathScripts}/compute.nb.elements.inbetween.pl --pathCoordinates1=${pathInteractions} --pathCoordinates2=${pathGenes} --type=genes --pathOutput=${pathResults}/statistics_${dataset}_interactions.txt
   
done

########################################################################
