#!/bin/bash

export user=${USER}

##########################################################################

if [ ${user} = "necsulea" ]; then
    export pathFinalData=/beegfs/data/necsulea/RegulatoryLandscapesManuscript
fi

export pathGO=${pathFinalData}/SupplementaryDataset3/gene_ontology

export release=94

##########################################################################

for sp in human mouse
do
    grep "multicellular organism development" ${pathGO}/${sp}_gene_ontology_Ensembl${release}.txt | cut -f 1 | sort -u >  ${pathGO}/${sp}_developmental_genes_Ensembl${release}.txt
done

##########################################################################
