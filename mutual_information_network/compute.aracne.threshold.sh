#!/bin/bash

export genome=$1
export cluster=$2

####################################################################################

if [ ${cluster} = "in2p3" ]; then
    export path=/sps/biometr/necsulea/RegulatoryLandscapes
    export pathTools=/sps/biometr/necsulea/Tools
fi

if [ ${cluster} = "cloud" ]; then
    export path=/mnt/RegulatoryLandscapes
    export pathTools=/mnt/Tools
fi

export pathResults=${path}/results/mutual_information_network/${genome}
export pathAracne=${pathTools}/ARACNe-AP/
export pathScripts=${path}/scripts/mutual_information_network

####################################################################################

export regulators=EnhancerList.txt

####################################################################################

java -Xmx20G -jar ${pathAracne}/dist/aracne.jar -e ${pathResults}/TPM.txt --tfs ${pathResults}/${regulators} -o ${pathResults}/${sp}/aracne_replicates/ --pvalue 1E-8 --seed 1 --calculateThreshold

####################################################################################
