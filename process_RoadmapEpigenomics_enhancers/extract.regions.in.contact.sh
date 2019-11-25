#!/bin/bash

###################################################################################

export pathAnouk=/beegfs/data/necsulea/RegulatoryLandscapes
export pathAlex=/beegfs/data/alaverre/Regulatory_landscape
export pathRoadmap=${pathAnouk}/data/RoadmapEpigenomics/hg38/regulatory_regions
export pathInteractions=${pathAlex}/result/Supplementary_dataset2_merged_fragments/human
export pathResults=${pathAnouk}/results/co_expression_analysis/RoadmapEpigenomics/
export pathScripts=${pathAnouk}/scripts/process_RoadmapEpigenomics_enhancers

###################################################################################

perl ${pathScripts}/extract.regions.in.contact.pl --pathContacts=${pathInteractions}/all_interactions_merged.txt --pathPromoters=${pathRoadmap}/promoter_regions/all_promoters.bed  --pathEnhancers=${pathRoadmap}/enhancer_regions/all_enhancers.bed --pathOutput=${pathResults}/promoters_enhancers_in_contact_real_data.txt

###################################################################################

perl ${pathScripts}/extract.regions.in.contact.pl --pathContacts=${pathInteractions}/simulations_merged.txt  --pathPromoters=${pathRoadmap}/promoter_regions/all_promoters.bed  --pathEnhancers=${pathRoadmap}/enhancer_regions/all_enhancers.bed --pathOutput=${pathResults}/promoters_enhancers_in_contact_simulated_data.txt

###################################################################################

