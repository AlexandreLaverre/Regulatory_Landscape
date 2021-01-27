#!/bin/bash

export sp=$1

########################################################################

export pathOriginal=/beegfs/data/necsulea/RegulatoryLandscapesManuscript/SupplementaryDataset1/${sp}
export pathAnnot=/beegfs/data/necsulea/RegulatoryLandscapesManuscript/SupplementaryDataset3/genes
export pathEnhancers=/beegfs/data/necsulea/RegulatoryLandscapesManuscript/SupplementaryDataset4/${sp}
export pathFragments=/beegfs/data/necsulea/RegulatoryLandscapesManuscript/SupplementaryDataset5/${sp}
export pathScripts=/beegfs/data/necsulea/RegulatoryLandscapes/scripts/sequence_composition

########################################################################

if [ ${sp} = "human" ]; then
    export pathBaits=${pathOriginal}/bait_coords_hg38.txt
fi

if [ ${sp} = "mouse" ]; then
    export pathBaits=${pathOriginal}/bait_coords_mm10.txt
fi

export pathGenes=${pathAnnot}/${sp}_genes_Ensembl94.txt

########################################################################

for dataset in original simulated
do
    ## fragments
    
    perl ${pathScripts}/compute.nb.elements.neighborhood.pl --pathCoordinates1=${pathFragments}/statistics_contacted_sequence_${dataset}.txt --pathCoordinates2=${pathBaits} --type=baits --pathOutput=${pathFragments}/statistics_contacted_sequence_${dataset}_with_nbbaits.txt
    
    perl ${pathScripts}/compute.nb.elements.neighborhood.pl --pathCoordinates1=${pathFragments}/statistics_contacted_sequence_${dataset}_with_nbbaits.txt --pathCoordinates2=${pathGenes} --type=genes --pathOutput=${pathFragments}/statistics_contacted_sequence_${dataset}_with_nbbaits_nbgenes.txt

    mv ${pathFragments}/statistics_contacted_sequence_${dataset}.txt ${pathFragments}/backup_statistics_contacted_sequence_${dataset}.txt
    mv ${pathFragments}/statistics_contacted_sequence_${dataset}_with_nbbaits_nbgenes.txt ${pathFragments}/statistics_contacted_sequence_${dataset}.txt

    rm ${pathFragments}/statistics_contacted_sequence_${dataset}_with_nbbaits.txt

    ## enhancers
    
    for enh in ENCODE FANTOM5 FOCS_GRO_seq RoadmapEpigenomics
    do
	if [ -e ${pathEnhancers}/${enh}/statistics_contacted_enhancers_${dataset}.txt ]; then
    	    perl ${pathScripts}/compute.nb.elements.neighborhood.pl --pathCoordinates1=${pathEnhancers}/${enh}/statistics_contacted_enhancers_${dataset}.txt  --pathCoordinates2=${pathBaits} --type=baits --pathOutput=${pathEnhancers}/${enh}/statistics_contacted_enhancers_${dataset}_with_nbbaits.txt

	    perl ${pathScripts}/compute.nb.elements.neighborhood.pl --pathCoordinates1=${pathEnhancers}/${enh}/statistics_contacted_enhancers_${dataset}_with_nbbaits.txt  --pathCoordinates2=${pathGenes} --type=genes --pathOutput=${pathEnhancers}/${enh}/statistics_contacted_enhancers_${dataset}_with_nbbaits_nbgenes.txt

	    
	    mv ${pathEnhancers}/${enh}/statistics_contacted_enhancers_${dataset}.txt ${pathEnhancers}/${enh}/backup_statistics_contacted_enhancers_${dataset}.txt
	    mv ${pathEnhancers}/${enh}/statistics_contacted_enhancers_${dataset}_with_nbbaits_nbgenes.txt ${pathEnhancers}/${enh}/statistics_contacted_enhancers_${dataset}.txt
	    
	    rm ${pathEnhancers}/${enh}/statistics_contacted_enhancers_${dataset}_with_nbbaits.txt

	fi
    done
done

########################################################################
