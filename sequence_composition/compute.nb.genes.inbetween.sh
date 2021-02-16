#!/bin/bash

export sp=$1

########################################################################

export pathSynteny=/beegfs/data/necsulea/RegulatoryLandscapesManuscript/SupplementaryDataset7/${sp}/synteny_conservation
export pathAnnot=/beegfs/data/necsulea/RegulatoryLandscapesManuscript/SupplementaryDataset3/genes
export pathScripts=/beegfs/data/necsulea/RegulatoryLandscapes/scripts/sequence_composition

########################################################################

export pathReferenceGenes=${pathAnnot}/${sp}_genes_Ensembl94.txt

########################################################################

for dataset in original simulated
do
    for enh in ENCODE  FANTOM5  FOCS_GRO_seq  RoadmapEpigenomics
    do
	for tg in human macaque mouse rat rabbit dog cow elephant opossum chicken
	do
	    if [ -e ${pathSynteny}/${enh}/${sp}2${tg}_${dataset}_synteny.txt ]; then
		
		if [ ${tg} = "macaque" ]; then
		    export pathTargetGenes=${pathAnnot}/${tg}_genes_Ensembl99.txt
		else
		    export pathTargetGenes=${pathAnnot}/${tg}_genes_Ensembl94.txt
		fi
		
		mv ${pathSynteny}/${enh}/${sp}2${tg}_${dataset}_synteny.txt ${pathSynteny}/${enh}/backup_${sp}2${tg}_${dataset}_synteny.txt

		perl ${pathScripts}/compute.nb.genes.inbetween.pl --pathContacts=${pathSynteny}/${enh}/backup_${sp}2${tg}_${dataset}_synteny.txt --pathReferenceGenes=${pathReferenceGenes} --pathTargetGenes=${pathTargetGenes} --pathOutput=${pathSynteny}/${enh}/${sp}2${tg}_${dataset}_synteny.txt

	    fi
	done
    done
done

########################################################################
