#!/bin/bash

export sp=$1
export cluster=$2

##############################################################

if [ ${cluster} = "pbil" ]; then
    export path=/bgeefs/data/necsulea/RegulatoryLandscapes
fi

if [ ${cluster} = "in2p3" ]; then
    export path=/sps/biometr/necsulea/RegulatoryLandscapes
fi

if [ ${cluster} = "cloud" ]; then
    export path=/mnt/RegulatoryLandscapes
fi

###########################################################################

export pathEnsembl=${path}/data/ensembl_annotations/${sp}
export pathStringTie=${path}/results/stringtie_assembly/${sp}/combined
export pathGenome=${path}/data/genome_sequences/${sp}
export pathScripts=${path}/scripts/expression_estimation

export release=94

##############################################################

for suffix in FilteredTranscripts_Ensembl${release} # AllTranscripts_Ensembl${release}
do
    perl ${pathScripts}/extract.cDNA.sequences.pl --pathAnnotGTF=${pathEnsembl}/${suffix}.gtf --forbiddenChromo=MT --pathGenomeSequence=${pathGenome}/genome_ensembl${release}.fa.gz --pathOutput=${pathEnsembl}/${suffix}_noMT.fa
done

###############################################################
