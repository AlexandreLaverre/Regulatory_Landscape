#!/bin/bash

export sp=$1

########################################################################

export pathSynteny=/beegfs/data/necsulea/RegulatoryLandscapesManuscript/SupplementaryDataset7/${sp}/synteny_conservation
export pathContacts=/beegfs/data/necsulea/RegulatoryLandscapesManuscript/SupplementaryDataset7/${sp}/contact_conservation
export pathFragmentInteractions=/beegfs/data/necsulea/RegulatoryLandscapesManuscript/SupplementaryDataset1/${sp}
export pathAnnot=/beegfs/data/necsulea/RegulatoryLandscapesManuscript/SupplementaryDataset3/transcripts
export pathScripts=/beegfs/data/necsulea/RegulatoryLandscapes/scripts/sequence_composition

########################################################################

export pathReferenceTranscriptCoords=${pathAnnot}/${sp}_transcripts_Ensembl94.txt

if [ ${sp} = "human" ]; then
    export pathBaitAnnotation=${pathFragmentInteractions}/bait_coords_hg38.txt 
fi

if [ ${sp} = "mouse" ]; then
    export pathBaitAnnotation=${pathFragmentInteractions}/bait_coords_mm10.txt 
fi

########################################################################

## synteny conservation

for dataset in original simulated
do
    for enh in ENCODE  FANTOM5  FOCS_GRO_seq  RoadmapEpigenomics
    do
	for tg in human macaque mouse rat rabbit dog cow elephant opossum chicken
	do
	    if [ ${tg} = "macaque" ]; then
		export pathTargetTranscriptCoords=${pathAnnot}/${tg}_transcripts_Ensembl99.txt
	    else
		export pathTargetTranscriptCoords=${pathAnnot}/${tg}_transcripts_Ensembl94.txt
	    fi
	    
	    if [ -e ${pathSynteny}/${enh}/${sp}2${tg}_${dataset}_synteny.txt ]; then
		mv ${pathSynteny}/${enh}/${sp}2${tg}_${dataset}_synteny.txt  ${pathSynteny}/${enh}/backup_${sp}2${tg}_${dataset}_synteny.txt 
		
		perl ${pathScripts}/compute.distance.TSS.enhancer.contacts.pl --pathContacts=${pathSynteny}/${enh}/backup_${sp}2${tg}_${dataset}_synteny.txt --pathReferenceBaitAnnotation=${pathBaitAnnotation} --pathReferenceTranscriptCoordinates=${pathReferenceTranscriptCoords}  --pathTargetTranscriptCoordinates=${pathTargetTranscriptCoords} --pathOutput=${pathSynteny}/${enh}/${sp}2${tg}_${dataset}_synteny.txt 
	    fi
	done
    done
done

########################################################################

## contact conservation

# for datasetref in original simulated
# do
#     for datasettg in original simulated
#     do
# 	for enh in ENCODE  FANTOM5  FOCS_GRO_seq RoadmapEpigenomics
# 	do
# 	    for tg in human mouse 
# 	    do	    
# 		export pathTargetTranscriptCoords=${pathAnnot}/${tg}_transcripts_Ensembl94.txt
		
# 		if [ -e ${pathContacts}/${enh}/${sp}_${datasetref}2${tg}_${datasettg}.txt ]; then
# 		    mv ${pathContacts}/${enh}/${sp}_${datasetref}2${tg}_${datasettg}.txt ${pathContacts}/${enh}/backup_${sp}_${datasetref}2${tg}_${datasettg}.txt
		    
# 		    perl ${pathScripts}/compute.distance.TSS.enhancer.contacts.pl --pathContacts=${pathContacts}/${enh}/backup_${sp}_${datasetref}2${tg}_${datasettg}.txt --pathReferenceBaitAnnotation=${pathBaitAnnotation} --pathReferenceTranscriptCoordinates=${pathReferenceTranscriptCoords} --pathTargetTranscriptCoordinates=${pathTargetTranscriptCoords} --pathOutput=${pathContacts}/${enh}/${sp}_${datasetref}2${tg}_${datasettg}.txt
# 		fi
# 	    done
# 	done
#     done
# done

########################################################################
