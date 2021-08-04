#!/bin/bash

export sp=$1
export dataset=$2
export cluster=$3

###################################################################################

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/necsulea/RegulatoryLandscapes
fi

export pathFOCS=${path}/data/FOCS/${sp}/${dataset}/
export pathResults=${path}/results/co_expression_analysis/${sp}/${dataset}
export pathScripts=${path}/scripts/promoter_enhancer_coexpression

###################################################################################

if [ -e ${pathResults}/expression_correlations_real_data.txt ]; then
    echo "already done"
else
    perl ${pathScripts}/compute.correlations.pl --pathContacts=${pathResults}/promoters_enhancers_in_contact_real_data.txt --pathPromoterExpression=${pathFOCS}/promoter_activity.txt --pathEnhancerExpression=${pathFOCS}/enhancer_activity.txt --pathOutput=${pathResults}/expression_correlations_real_data.txt
fi

###################################################################################

if [ -e ${pathResults}/expression_correlations_simulated_data.txt ]; then
    echo "already done"
else
    perl ${pathScripts}/compute.correlations.pl --pathContacts=${pathResults}/promoters_enhancers_in_contact_simulated_data.txt --pathPromoterExpression=${pathFOCS}/promoter_activity.txt --pathEnhancerExpression=${pathFOCS}/enhancer_activity.txt --pathOutput=${pathResults}/expression_correlations_simulated_data.txt
fi

###################################################################################

if [ -e ${pathResults}/expression_correlations_neigbor_promoters_enhancers.txt ]; then
    echo "already done"
else
    perl ${pathScripts}/compute.correlations.pl --pathContacts=${pathResults}/neighbor_promoters_enhancers.txt --pathPromoterExpression=${pathFOCS}/promoter_activity.txt --pathEnhancerExpression=${pathFOCS}/enhancer_activity.txt --pathOutput=${pathResults}/expression_correlations_neighbor_promoters_enhancers.txt
fi

###################################################################################
