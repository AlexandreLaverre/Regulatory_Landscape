#!/bin/bash

###################################################################################

export path=/beegfs/data/necsulea/RegulatoryLandscapes
export pathFOCS=${path}/data/FOCS
export pathCoExpression=${path}/results/co_expression_analysis/
export pathManuscript=/beegfs/data/necsulea/RegulatoryLandscapesManuscript
export pathResults=${pathManuscript}/SupplementaryDataset8

###################################################################################

for sp in human mouse
do
    for dataset in ENCODE  FANTOM5  GRO-seq  RoadmapEpigenomics
    do
	if [ ${dataset} = "GRO-seq" ]; then
	    outdataset=FOCS_GRO_seq
	else
	    outdataset=${dataset}
	fi

	if [ -e ${pathManuscript}/${sp}/${outdataset} ]; then
	    cp ${pathFOCS}/${sp}/${dataset}/promoter_coords_activity.txt ${pathManuscript}/${sp}/${outdataset}/promoter_activity.txt
	    cp ${pathFOCS}/${sp}/${dataset}/enhancer_activity.txt ${pathManuscript}/${sp}/${outdataset}/enhancer_activity.txt

	    cp ${pathCoExpression}/${sp}/${dataset}/*  ${pathManuscript}/${sp}/${outdataset}/
	fi
    done
done

###################################################################################
