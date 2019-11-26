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

if [ -e ${pathResults}/all_expression_correlations_promoters_enhancers_in_contact_real_data.txt ]; then
    echo "already done"
else
    perl ${pathScripts}/compute.correlations.pl --pathContacts=${pathResults}/all_promoters_enhancers_in_contact_real_data.txt --pathPromoterExpression=${pathRoadmapHg19}/promoter_regions/all_promoters_coverage_allsamples.txt --pathEnhancerExpression=${pathRoadmapHg19}/enhancer_regions/all_enhancers_coverage_allsamples.txt --pathOutput=${pathResults}/all_expression_correlations_promoters_enhancers_in_contact_real_data.txt
fi

if [ -e ${pathResults}/selected_expression_correlations_promoters_enhancers_in_contact_real_data.txt ]; then
    echo "already done"
else
    perl ${pathScripts}/compute.correlations.pl --pathContacts=${pathResults}/selected_promoters_enhancers_in_contact_real_data.txt --pathPromoterExpression=${pathRoadmapHg19}/promoter_regions/all_promoters_coverage_allsamples.txt --pathEnhancerExpression=${pathRoadmapHg19}/enhancer_regions/all_enhancers_coverage_allsamples.txt --pathOutput=${pathResults}/selected_expression_correlations_promoters_enhancers_in_contact_real_data.txt
fi

if [ -e ${pathResults}/FOCS_expression_correlations_promoters_enhancers_in_contact_real_data.txt ]; then
    echo "already done"
else
    perl ${pathScripts}/compute.correlations.pl --pathContacts=${pathResults}/FOCS_promoters_enhancers_in_contact_real_data.txt --pathPromoterExpression=${pathRoadmapHg19}/promoter_regions/promoters_FOCS_coverage_allsamples.txt --pathEnhancerExpression=${pathRoadmapHg19}/enhancer_regions/enhancers_FOCS_coverage_allsamples.txt --pathOutput=${pathResults}/FOCS_expression_correlations_promoters_enhancers_in_contact_real_data.txt
fi

###################################################################################

if [ -e ${pathResults}/all_expression_correlations_promoters_enhancers_in_contact_simulated_data.txt ]; then
    echo "already done"
else
    perl ${pathScripts}/compute.correlations.pl --pathContacts=${pathResults}/all_promoters_enhancers_in_contact_simulated_data.txt --pathPromoterExpression=${pathRoadmapHg19}/promoter_regions/all_promoters_coverage_allsamples.txt --pathEnhancerExpression=${pathRoadmapHg19}/enhancer_regions/all_enhancers_coverage_allsamples.txt --pathOutput=${pathResults}/all_expression_correlations_promoters_enhancers_in_contact_simulated_data.txt
fi

if [ -e ${pathResults}/selected_expression_correlations_promoters_enhancers_in_contact_simulated_data.txt ]; then
    echo "already done"
else
    perl ${pathScripts}/compute.correlations.pl --pathContacts=${pathResults}/selected_promoters_enhancers_in_contact_simulated_data.txt --pathPromoterExpression=${pathRoadmapHg19}/promoter_regions/all_promoters_coverage_allsamples.txt --pathEnhancerExpression=${pathRoadmapHg19}/enhancer_regions/all_enhancers_coverage_allsamples.txt --pathOutput=${pathResults}/selected_expression_correlations_promoters_enhancers_in_contact_simulated_data.txt
fi


if [ -e ${pathResults}/FOCS_expression_correlations_promoters_enhancers_in_contact_simulated_data.txt ]; then
    echo "already done"
else
    perl ${pathScripts}/compute.correlations.pl --pathContacts=${pathResults}/FOCS_promoters_enhancers_in_contact_simulated_data.txt --pathPromoterExpression=${pathRoadmapHg19}/promoter_regions/promoters_FOCS_coverage_allsamples.txt --pathEnhancerExpression=${pathRoadmapHg19}/enhancer_regions/enhancers_FOCS_coverage_allsamples.txt --pathOutput=${pathResults}/FOCS_expression_correlations_promoters_enhancers_in_contact_simulated_data.txt
fi

###################################################################################
