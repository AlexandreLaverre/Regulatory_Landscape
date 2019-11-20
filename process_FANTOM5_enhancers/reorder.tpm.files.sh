#!/bin/bash

export cluster=$1

##################################################################

if [ ${cluster} = "in2p3" ]; then
    export path=/sps/biometr/necsulea/RegulatoryLandscapes
fi

export pathData=${path}/data/FANTOM5
export pathScripts=${path}/scripts/process_FANTOM5_enhancers 

##################################################################

for genome in hg19 mm9
do
    if [ ${genome} = "hg19" ]; then
	export sp="human"
    fi

    if [ ${genome} = "mm9" ]; then
	export sp="mouse"
    fi
    
    perl ${pathScripts}/reorder.tpm.files.pl --pathSelectedPeaks=${pathData}/${genome}/SelectedCAGEPeaks.txt --pathSelectedSamples=${pathData}/${genome}/SelectedSamples.txt  --pathOriginalTPM=${pathData}/${genome}/${genome}.cage_peak_phase1and2combined_tpm.osc.txt.gz  --pathReorderedTPM=${pathData}/${genome}/${genome}.cage_peak_phase1and2combined_tpm.osc.reordered.txt

    perl ${pathScripts}/reorder.tpm.files.pl --pathSelectedPeaks=${pathData}/${genome}/SelectedGenes.txt --pathSelectedSamples=${pathData}/${genome}/SelectedSamples.txt  --pathOriginalTPM=${pathData}/${genome}/${genome}.gene_phase1and2combined_tpm.osc.txt.gz   --pathReorderedTPM=${pathData}/${genome}/${genome}.gene_phase1and2combined_tpm.osc.reordered.txt
    
    perl ${pathScripts}/reorder.tpm.files.pl --pathSelectedPeaks=${pathData}/${genome}/SelectedEnhancers.txt --pathSelectedSamples=${pathData}/${genome}/SelectedSamples.txt --pathOriginalTPM=${pathData}/${genome}/${sp}_permissive_enhancers_phase_1_and_2_expression_tpm_matrix.txt.gz --pathReorderedTPM=${pathData}/${genome}/${sp}_permissive_enhancers_phase_1_and_2_expression_tpm_matrix.reordered.txt 
done

##################################################################
