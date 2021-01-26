#!/bin/bash

export sp=$1

########################################################################

export pathGenomes=/beegfs/data/alaverre/Regulatory_landscape/data/genome_ref
export pathEnhancers=/beegfs/data/necsulea/RegulatoryLandscapesManuscript/SupplementaryDataset4/${sp}
export pathFragments=/beegfs/data/necsulea/RegulatoryLandscapesManuscript/SupplementaryDataset5/${sp}
export pathScripts=/beegfs/data/necsulea/RegulatoryLandscapes/scripts/sequence_composition

########################################################################

if [ ${sp} = "human" ]; then
    export pathGenomeSequence=${pathGenomes}/Homo_sapiens.GRCh38.dna.primary_assembly.fa 
fi

if [ ${sp} = "mouse" ]; then
    export pathGenomeSequence=${pathGenomes}/Mus_musculus.mm10.dna.primary_assembly.fa
fi

########################################################################

## fragments

perl ${pathScripts}/compute.GC.content.pl --pathGenomeSequence=${pathGenomeSequence} --pathCoordinates=${pathFragments}/statistics_contacted_sequence_original.txt --pathOutput=${pathFragments}/statistics_contacted_sequence_original_withGCcontent.txt

perl ${pathScripts}/compute.GC.content.pl --pathGenomeSequence=${pathGenomeSequence} --pathCoordinates=${pathFragments}/statistics_contacted_sequence_simulated.txt --pathOutput=${pathFragments}/statistics_contacted_sequence_simulated_withGCcontent.txt

########################################################################

## enhancers

for enh in ENCODE FANTOM5 FOCS_GRO_seq RoadmapEpigenomics
do
    if [ -e ${pathEnhancers}/${enh}/statistics_contacted_enhancers_original.txt ]; then
	perl ${pathScripts}/compute.GC.content.pl --pathGenomeSequence=${pathGenomeSequence} --pathCoordinates=${pathEnhancers}/${enh}/statistics_contacted_enhancers_original.txt --pathOutput=${pathEnhancers}/${enh}/statistics_contacted_enhancers_original_withGCcontent.txt

	perl ${pathScripts}/compute.GC.content.pl --pathGenomeSequence=${pathGenomeSequence} --pathCoordinates=${pathEnhancers}/${enh}/statistics_contacted_enhancers_simulated.txt --pathOutput=${pathEnhancers}/${enh}/statistics_contacted_enhancers_simulated_withGCcontent.txt
    fi
done

########################################################################
