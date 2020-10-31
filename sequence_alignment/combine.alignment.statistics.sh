#!/bin/bash

export ref=$1
export type=$2
export enh=$3

#####################################################################

export pathFinalData=/beegfs/data/necsulea/RegulatoryLandscapesManuscript
export pathInteractions=${pathFinalData}/SupplementaryDataset1/${ref}
export pathEnhancers=${pathFinalData}/SupplementaryDataset4/${ref}
export pathSequenceConservation=${pathFinalData}/SupplementaryDataset7/${ref}/sequence_conservation
export pathScripts=/beegfs/data/necsulea/RegulatoryLandscapes/scripts/sequence_alignment

#####################################################################

export tgSpecies=""
export pathsAln=""

if [ ${ref} = "human" ]; then
    export genome="hg38"
    
    for sp in macaque mouse rat rabbit cow dog elephant opossum chicken
    do
	export tgSpecies=${sp},${tgSpecies}

	if [ ${type} = "enhancers" ]; then
	    export pathsAln=${pathSequenceConservation}/${type}/${enh}/AlignmentStatistics_Excluding_Exons_${ref}2${sp}.txt,${pathsAln}
	fi

	if [ ${type} = "restriction_fragments" ]; then
	    export pathsAln=${pathSequenceConservation}/${type}/AlignmentStatistics_Excluding_Exons_${ref}2${sp}.txt,${pathsAln}
	fi

    done
fi

if [ ${ref} = "mouse" ]; then
    export genome="mm10"
    
    for sp in macaque human rat rabbit cow dog elephant opossum chicken
    do
	export tgSpecies=${sp},${tgSpecies}
	
	if [ ${type} = "enhancers" ]; then
	    export pathsAln=${pathSequenceConservation}/${type}/${enh}/AlignmentStatistics_Excluding_Exons_${ref}2${sp}.txt,${pathsAln}
	fi

	if [ ${type} = "restriction_fragments" ]; then
	    export pathsAln=${pathSequenceConservation}/${type}/AlignmentStatistics_Excluding_Exons_${ref}2${sp}.txt,${pathsAln}
	fi
    done
fi

#####################################################################

if [ ${type} = "enhancers" ]; then
    export pathCoords=${pathEnhancers}/${enh}/enhancer_coordinates.bed
    export pathResults=${pathSequenceConservation}/${type}/${enh}
fi

if [ ${type} = "restriction_fragments" ]; then
    export pathCoords=${pathInteractions}/frag_coords_${genome}.txt
    export pathResults=${pathSequenceConservation}/${type}
fi

#####################################################################

perl ${pathScripts}/combine.alignment.statistics.pl --pathCoordinates=${pathCoords} --refSpecies=${ref} --targetSpecies=${tgSpecies} --pathsAlignmentStatistics=${pathsAln} --pathOutputUngappedAlignment=${pathResults}/AlignmentStatistics_Excluding_Exons_UngappedAlignment_AllSpecies.txt --pathOutputSequenceIdentity=${pathResults}/AlignmentStatistics_Excluding_Exons_IdenticalSequence_AllSpecies.txt  

#####################################################################
