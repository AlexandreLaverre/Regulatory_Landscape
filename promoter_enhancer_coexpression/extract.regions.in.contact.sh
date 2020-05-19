#!/bin/bash

export sp=$1
export dataset=$2

###################################################################################

export pathAnouk=/beegfs/data/necsulea/RegulatoryLandscapes
export pathAlex=/beegfs/data/alaverre/Regulatory_landscape
export pathFOCS=${pathAnouk}/data/FOCS/${sp}/${dataset}/
export pathOriginalInteractions=${pathAlex}/result/Supplementary_dataset1_original_interactions/${sp}
export pathSimulatedInteractions=${pathAlex}/result/Supplementary_dataset2_simulated_interactions/${sp}
export pathResults=${pathAnouk}/results/co_expression_analysis/${sp}/${dataset}
export pathScripts=${pathAnouk}/scripts/promoter_enhancer_coexpression

###################################################################################

if [ -e ${pathResults} ]; then
    echo "output path already there"
else
    mkdir -p ${pathResults} 
fi

###################################################################################

if [ ${sp} = "human" ]; then
    export assembly="hg38"
fi

if [ ${sp} = "mouse" ]; then
    export assembly="mm10"
fi

###################################################################################

if [ -e ${pathResults}/promoters_enhancers_in_contact_real_data.txt ]; then
    echo "already done"
else
    perl ${pathScripts}/extract.regions.in.contact.pl --pathContacts=${pathOriginalInteractions}/all_interactions.txt --pathPromoters=${pathFOCS}/promoter_coordinates_${assembly}.bed  --pathEnhancers=${pathFOCS}/enhancer_coordinates_${assembly}.bed --pathOutput=${pathResults}/promoters_enhancers_in_contact_real_data.txt
fi

###################################################################################

if [ -e ${pathResults}/promoters_enhancers_in_contact_simulated_data.txt ]; then
    echo "already done"
else
    perl ${pathScripts}/extract.regions.in.contact.pl --pathContacts=${pathSimulatedInteractions}/simulated_all_interactions.txt --pathPromoters=${pathFOCS}/promoter_coordinates_${assembly}.bed  --pathEnhancers=${pathFOCS}/enhancer_coordinates_${assembly}.bed --pathOutput=${pathResults}/promoters_enhancers_in_contact_simulated_data.txt
fi

###################################################################################
