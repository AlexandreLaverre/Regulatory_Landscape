#!/bin/bash

export sp=$1
export dataset=$2

###################################################################################

export pathAnouk=/beegfs/data/necsulea/RegulatoryLandscapes
export pathAlex=/beegfs/data/alaverre/Regulatory_landscape
export pathFOCS=${pathAnouk}/data/FOCS/${sp}/${dataset}/
export pathResults=${pathAnouk}/results/co_expression_analysis/${sp}/${dataset}
export pathScripts=${pathAnouk}/scripts/promoter_enhancer_coexpression

###################################################################################

if [ -e ${pathResults}/expression_correlations_real_data.txt ]; then
    echo "already done"
else
    perl ${pathScripts}/compute.correlations.pl --pathContacts=${pathResults}/promoters_enhancers_in_contact_real_data.txt --pathPromoterExpression=${pathFOCS}/promoter_activity.txt --pathEnhancerExpression=${pathFOCS}/enhancer_activity.txt --pathOutput=${pathResults}/expression_correlations_real_data.txt
fi

###################################################################################
