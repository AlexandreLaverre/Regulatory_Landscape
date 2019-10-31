#!/bin/bash

export genome=$1
export cloud=$2

##################################################################

if [ ${cloud} = "in2p3" ]; then
    export path=/sps/biometr/necsulea/RegulatoryLandscapes
fi

export pathData=${path}/data/FANTOM5/${genome}
export pathResults=${path}/results/mutual_information_network/${genome}
export pathScripts=${path}/scripts/process_FANTOM5_enhancers 

##################################################################

if [ ${genome} = "hg19" ]; then
    export sp="human"
fi

if [ ${genome} = "mm9" ]; then
    export sp="mouse"
fi

##################################################################

cp ${pathData}/${sp}_permissive_enhancers_phase_1_and_2_expression_tpm_matrix.reordered.txt ${pathResults}/TPM.txt
sed '1,3d' ${pathData}/${genome}.gene_phase1and2combined_tpm.osc.reordered.txt >> ${pathResults}/TPM.txt

cut -f 1 ${pathData}/${sp}_permissive_enhancers_phase_1_and_2_expression_tpm_matrix.reordered.txt | sed '1d' > ${pathResults}/Enhancers.txt

##################################################################
