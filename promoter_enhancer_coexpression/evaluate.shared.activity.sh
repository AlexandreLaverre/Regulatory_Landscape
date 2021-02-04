#!/bin/bash

export sp=$1
export dataset=$2
export minTPM=$3

###################################################################################

export path=/beegfs/data/necsulea/RegulatoryLandscapes
export pathFOCS=${path}/data/FOCS/${sp}/${dataset}/
export pathResults=${path}/results/co_expression_analysis/${sp}/${dataset}
export pathScripts=${path}/scripts/promoter_enhancer_coexpression

###################################################################################

for data in real simulated
do
    if [ -e ${pathResults}/shared_activity_minTPM${minTPM}_${data}_data.txt ]; then
	echo "already done"
    else
	perl ${pathScripts}/evaluate.shared.activity.pl --pathContacts=${pathResults}/promoters_enhancers_in_contact_${data}_data.txt --pathPromoterExpression=${pathFOCS}/promoter_activity.txt --pathEnhancerExpression=${pathFOCS}/enhancer_activity.txt --minValue=${minTPM} --pathOutput=${pathResults}/shared_activity_minTPM${minTPM}_${data}_data.txt
    fi
done

###################################################################################
