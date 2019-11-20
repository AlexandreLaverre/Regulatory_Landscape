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
    
    perl ${pathScripts}/compute.average.tpm.pl --pathSampleCategories=${pathData}/${genome}/SampleTypes.txt --pathOriginalTPM=${pathData}/${genome}/${genome}.gene_phase1and2combined_tpm.osc.txt.gz  --pathAverageTPM=${pathData}/${genome}/${genome}.gene_phase1and2combined_tpm.osc.averagevalues.txt
    
    perl ${pathScripts}/compute.average.tpm.pl --pathSampleCategories=${pathData}/${genome}/SampleTypes.txt --pathOriginalTPM=${pathData}/${genome}/${sp}_permissive_enhancers_phase_1_and_2_expression_tpm_matrix.txt.gz --pathAverageTPM=${pathData}/${genome}/${sp}_permissive_enhancers_phase_1_and_2_expression_tpm_matrix.averagevalues.txt 
done

##################################################################
