#!/bin/bash

###################################################################################

export pathAnouk=/beegfs/data/necsulea/RegulatoryLandscapes
export pathAlex=/beegfs/data/alaverre/Regulatory_landscape
export pathRoadmap=${pathAnouk}/data/RoadmapEpigenomics/hg38/regulatory_regions
export pathRoadmapHg19=${pathAnouk}/data/RoadmapEpigenomics/hg19/regulatory_regions
export pathInteractions=${pathAlex}/result/Supplementary_dataset2_merged_fragments/human
export pathResults=${pathAnouk}/results/co_expression_analysis/RoadmapEpigenomics/
export pathScripts=${pathAnouk}/scripts/process_RoadmapEpigenomics_enhancers

###################################################################################

if [ -e ${pathResults}/expression_correlations_promoters_enhancers_in_contact_real_data.txt ]; then
    echo "already done"
else
    perl ${pathScripts}/compute.correlations.pl --pathContacts=${pathResults}/promoters_enhancers_in_contact_real_data.txt --pathPromoterExpression=${pathRoadmapHg19}/promoter_regions/all_promoters_coverage_allsamples.txt --pathEnhancerExpression=${pathRoadmapHg19}/enhancer_regions/all_enhancers_coverage_allsamples.txt --pathOutput=${pathResults}/expression_correlations_promoters_enhancers_in_contact_real_data.txt
fi

if [ -e ${pathResults}/expression_correlations_promoters_enhancers_in_contact_real_data_selected_regions.txt ]; then
    echo "already done"
else
    perl ${pathScripts}/compute.correlations.pl --pathContacts=${pathResults}/promoters_enhancers_in_contact_real_data_selected_regions.txt --pathPromoterExpression=${pathRoadmapHg19}/promoter_regions/all_promoters_coverage_allsamples.txt --pathEnhancerExpression=${pathRoadmapHg19}/enhancer_regions/all_enhancers_coverage_allsamples.txt --pathOutput=${pathResults}/expression_correlations_promoters_enhancers_in_contact_real_data_selected_regions.txt
fi


###################################################################################

if [ -e ${pathResults}/expression_correlations_promoters_enhancers_in_contact_simulated_data.txt ]; then
    echo "already done"
else
    perl ${pathScripts}/compute.correlations.pl --pathContacts=${pathResults}/promoters_enhancers_in_contact_simulated_data.txt --pathPromoterExpression=${pathRoadmapHg19}/promoter_regions/all_promoters_coverage_allsamples.txt --pathEnhancerExpression=${pathRoadmapHg19}/enhancer_regions/all_enhancers_coverage_allsamples.txt --pathOutput=${pathResults}/expression_correlations_promoters_enhancers_in_contact_simulated_data.txt
fi



if [ -e ${pathResults}/expression_correlations_promoters_enhancers_in_contact_simulated_data_selected_regions.txt ]; then
    echo "already done"
else
    perl ${pathScripts}/compute.correlations.pl --pathContacts=${pathResults}/promoters_enhancers_in_contact_simulated_data_selected_regions.txt --pathPromoterExpression=${pathRoadmapHg19}/promoter_regions/all_promoters_coverage_allsamples.txt --pathEnhancerExpression=${pathRoadmapHg19}/enhancer_regions/all_enhancers_coverage_allsamples.txt --pathOutput=${pathResults}/expression_correlations_promoters_enhancers_in_contact_simulated_data_selected_regions.txt
fi


###################################################################################
