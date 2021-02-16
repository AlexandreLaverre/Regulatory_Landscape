#!/bin/bash

export sp=$1

########################################################################

export pathContacts=/beegfs/data/necsulea/RegulatoryLandscapesManuscript/SupplementaryDataset4/${sp}
export pathFragmentInteractions=/beegfs/data/necsulea/RegulatoryLandscapesManuscript/SupplementaryDataset1/${sp}
export pathAnnot=/beegfs/data/necsulea/RegulatoryLandscapesManuscript/SupplementaryDataset3/transcripts
export pathScripts=/beegfs/data/necsulea/RegulatoryLandscapes/scripts/sequence_composition

########################################################################

export pathTranscriptCoords=${pathAnnot}/${sp}_transcripts_Ensembl94.txt

if [ ${sp} = "human" ]; then
    export pathBaitAnnotation=${pathFragmentInteractions}/bait_coords_hg38.txt 
fi

if [ ${sp} = "mouse" ]; then
    export pathBaitAnnotation=${pathFragmentInteractions}/bait_coords_mm10.txt 
fi

########################################################################

for dataset in original simulated
do
    for enh in ENCODE  FANTOM5  FOCS_GRO_seq  RoadmapEpigenomics
    do

	if [ -e ${pathContacts}/${enh}/backup_gene_enhancer_contacts_${dataset}_interactions.txt ]; then
	    echo "already backed-up, not doing anything"
	else
	    if [ -e ${pathContacts}/${enh}/gene_enhancer_contacts_${dataset}_interactions.txt ]; then
		mv ${pathContacts}/${enh}/gene_enhancer_contacts_${dataset}_interactions.txt ${pathContacts}/${enh}/backup_gene_enhancer_contacts_${dataset}_interactions.txt
		
		perl ${pathScripts}/compute.distance.TSS.enhancer.pl --pathContacts=${pathContacts}/${enh}/backup_gene_enhancer_contacts_${dataset}_interactions.txt --pathBaitAnnotation=${pathBaitAnnotation} --pathTranscriptCoordinates=${pathTranscriptCoords} --pathOutput=${pathContacts}/${enh}/gene_enhancer_contacts_${dataset}_interactions.txt
	    fi
	fi
    done
done

########################################################################
