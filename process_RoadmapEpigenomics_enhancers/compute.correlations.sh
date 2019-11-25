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

perl ${pathScripts}/compute.correlations.pl --pathContacts=${pathResults}/promoters_enhancers_in_contact_real_data.txt --pathPromoterExpression=${pathRoadmapHg19}/promoter_regions/all_promoters_coverage_allsamples.txt --pathEnhancerExpression=${pathRoadmapHg19}/enhancer_regions/all_enhancers_coverage_allsamples.txt --pathOutput=${pathResults}/expression_correlations_promoters_enhancers_in_contact_real_data.txt

###################################################################################

perl ${pathScripts}/compute.correlations.pl --pathContacts=${pathResults}/promoters_enhancers_in_contact_simulated_data.txt --pathPromoterExpression=${pathRoadmapHg19}/promoter_regions/all_promoters_coverage_allsamples.txt --pathEnhancerExpression=${pathRoadmapHg19}/enhancer_regions/all_enhancers_coverage_allsamples.txt --pathOutput=${pathResults}/expression_correlations_promoters_enhancers_in_contact_simulated_data.txt

###################################################################################
