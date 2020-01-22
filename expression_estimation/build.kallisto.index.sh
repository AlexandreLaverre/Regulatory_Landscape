#!/bin/bash

export species=$1
export type=$2
export cluster=$3

####################################################################################

if [ ${cluster} = "pbil" ]; then
    export path=/bgeefs/data/necsulea/RegulatoryLandscapes
fi

if [ ${cluster} = "in2p3" ]; then
    export path=/sps/biometr/necsulea/RegulatoryLandscapes
fi

if [ ${cluster} = "cloud" ]; then
    export path=/mnt/RegulatoryLandscapes
fi

####################################################################################

export pathResults=${path}/results/kallisto_indexes/${species}
export pathAnnot=${path}/data/ensembl_annotations/${species}
export pathScripts=${path}/scripts/expression_estimation

export release=94

####################################################################################


if [ ${type} = "FilteredTranscripts" ]; then
    export pathFasta=${pathAnnot}/FilteredTranscripts_Ensembl${release}_noMT.fa
fi


if [ ${type} = "AllTranscripts" ]; then
    export pathFasta=${pathAnnot}/AllTranscripts_Ensembl${release}_noMT.fa
fi

####################################################################################

if [ -e ${pathResults} ]; then
    echo "path results already there"
else
    mkdir ${pathResults}
fi

####################################################################################

kallisto index -i ${pathResults}/${type} ${pathFasta}

####################################################################################
