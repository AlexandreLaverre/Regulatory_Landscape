#!/bin/bash

export sp=$1
export dataset=$2

#######################################################################

export path=/beegfs/data/necsulea/RegulatoryLandscapes
export pathLiftOver=${path}/data/liftOver_files
export pathFOCS=${path}/data/FOCS

#######################################################################

if [ ${sp} = "human" ]; then
    export oldversion=hg19
    export newversion=hg38
    export alignmentFile=${pathLiftOver}/hg19ToHg38.over.chain.gz
fi

if [ ${sp} = "mouse" ]; then
    export oldversion=mm9
    export newversion=mm10
    export alignmentFile=${pathLiftOver}/mm9ToMm10.over.chain.gz
fi

#######################################################################

liftOver ${pathFOCS}/${sp}/${dataset}/promoter_coordinates_${oldversion}.bed ${alignmentFile} ${pathFOCS}/${sp}/${dataset}/promoter_coordinates_${newversion}.bed ${pathFOCS}/${sp}/${dataset}/promoter_coordinates_${newversion}.unmapped

liftOver ${pathFOCS}/${sp}/${dataset}/enhancer_coordinates_${oldversion}.bed ${alignmentFile} ${pathFOCS}/${sp}/${dataset}/enhancer_coordinates_${newversion}.bed ${pathFOCS}/${sp}/${dataset}/enhancer_coordinates_${newversion}.unmapped

#######################################################################
