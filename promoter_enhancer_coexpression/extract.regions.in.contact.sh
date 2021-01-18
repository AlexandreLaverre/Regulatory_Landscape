#!/bin/bash

export sp=$1
export dataset=$2

###################################################################################

export path=/beegfs/data/necsulea/RegulatoryLandscapes
export pathManuscript=/beegfs/data/necsulea/RegulatoryLandscapesManuscript
export pathFOCS=${path}/data/FOCS/${sp}/${dataset}/
export pathOriginalInteractions=${pathManuscript}/SupplementaryDataset1/${sp}
export pathSimulatedInteractions=${pathManuscript}/SupplementaryDataset2/${sp}
export pathResults=${path}/results/co_expression_analysis/${sp}/${dataset}
export pathScripts=${path}/scripts/promoter_enhancer_coexpression

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
