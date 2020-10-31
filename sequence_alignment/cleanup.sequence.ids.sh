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

for sp in human macaque mouse rat rabbit cow dog elephant opossum chicken
do
    if [ ${type} = "enhancers" ]; then
	if [ -e ${pathSequenceConservation}/${type}/${enh}/AlignmentStatistics_Excluding_Exons_${ref}2${sp}.txt ]; then
	    mv ${pathSequenceConservation}/${type}/${enh}/AlignmentStatistics_Excluding_Exons_${ref}2${sp}.txt ${pathSequenceConservation}/${type}/${enh}/AlignmentStatistics_Excluding_Exons_${ref}2${sp}_backup.txt
	fi

	if [ -e ${pathSequenceConservation}/${type}/${enh}/AlignmentStatistics_Excluding_Exons_${ref}2${sp}_backup.txt ]; then
	    perl ${pathScripts}/cleanup.sequence.ids.pl --pathInput=${pathSequenceConservation}/${type}/${enh}/AlignmentStatistics_Excluding_Exons_${ref}2${sp}_backup.txt --pathOutput=${pathSequenceConservation}/${type}/${enh}/AlignmentStatistics_Excluding_Exons_${ref}2${sp}.txt
	fi

    fi

    if [ ${type} = "restriction_fragments" ]; then
	if [ -e ${pathSequenceConservation}/${type}/AlignmentStatistics_Excluding_Exons_${ref}2${sp}.txt ]; then
	    mv ${pathSequenceConservation}/${type}/AlignmentStatistics_Excluding_Exons_${ref}2${sp}.txt ${pathSequenceConservation}/${type}/AlignmentStatistics_Excluding_Exons_${ref}2${sp}_backup.txt
	fi

	if [ -e ${pathSequenceConservation}/${type}/AlignmentStatistics_Excluding_Exons_${ref}2${sp}_backup.txt ]; then
	    perl ${pathScripts}/cleanup.sequence.ids.pl --pathInput=${pathSequenceConservation}/${type}/AlignmentStatistics_Excluding_Exons_${ref}2${sp}_backup.txt --pathOutput=${pathSequenceConservation}/${type}/AlignmentStatistics_Excluding_Exons_${ref}2${sp}.txt
	fi
    fi
    
done


#####################################################################
#####################################################################
