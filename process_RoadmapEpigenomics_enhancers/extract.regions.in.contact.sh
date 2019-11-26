#!/bin/bash

###################################################################################

export pathAnouk=/beegfs/data/necsulea/RegulatoryLandscapes
export pathAlex=/beegfs/data/alaverre/Regulatory_landscape
export pathRoadmap=${pathAnouk}/data/RoadmapEpigenomics/hg38/regulatory_regions
export pathInteractions=${pathAlex}/result/Supplementary_dataset2_merged_fragments/human
export pathResults=${pathAnouk}/results/co_expression_analysis/RoadmapEpigenomics/
export pathScripts=${pathAnouk}/scripts/process_RoadmapEpigenomics_enhancers

###################################################################################

if [ -e ${pathResults}/all_promoters_enhancers_in_contact_real_data.txt ]; then
    echo "already done"
else
    perl ${pathScripts}/extract.regions.in.contact.pl --pathContacts=${pathInteractions}/all_interactions_merged.txt --pathPromoters=${pathRoadmap}/promoter_regions/all_promoters.bed  --pathEnhancers=${pathRoadmap}/enhancer_regions/all_enhancers.bed --pathOutput=${pathResults}/all_promoters_enhancers_in_contact_real_data.txt
fi


if [ -e ${pathResults}/selected_promoters_enhancers_in_contact_real_data.txt ]; then
    echo "already done"
else
    perl ${pathScripts}/extract.regions.in.contact.pl --pathContacts=${pathInteractions}/all_interactions_merged.txt --pathPromoters=${pathRoadmap}/promoter_regions/selected_promoters.bed  --pathEnhancers=${pathRoadmap}/enhancer_regions/selected_enhancers.bed --pathOutput=${pathResults}/selected_promoters_enhancers_in_contact_real_data.txt
fi


if [ -e ${pathResults}/FOCS_promoters_enhancers_in_contact_real_data.txt ]; then
    echo "already done"
else
    perl ${pathScripts}/extract.regions.in.contact.pl --pathContacts=${pathInteractions}/all_interactions_merged.txt --pathPromoters=${pathRoadmap}/promoter_regions/promoters_FOCS.bed  --pathEnhancers=${pathRoadmap}/enhancer_regions/enhancers_FOCS.bed --pathOutput=${pathResults}/FOCS_promoters_enhancers_in_contact_real_data.txt
fi

###################################################################################

if [ -e ${pathResults}/all_promoters_enhancers_in_contact_simulated_data.txt ]; then
    echo "already done"
else
    perl ${pathScripts}/extract.regions.in.contact.pl --pathContacts=${pathInteractions}/simulations_merged.txt  --pathPromoters=${pathRoadmap}/promoter_regions/all_promoters.bed  --pathEnhancers=${pathRoadmap}/enhancer_regions/all_enhancers.bed --pathOutput=${pathResults}/all_promoters_enhancers_in_contact_simulated_data.txt
fi


if [ -e ${pathResults}/selected_promoters_enhancers_in_contact_simulated_data.txt ]; then
    echo "already done"
else
    perl ${pathScripts}/extract.regions.in.contact.pl --pathContacts=${pathInteractions}/simulations_merged.txt  --pathPromoters=${pathRoadmap}/promoter_regions/selected_promoters.bed  --pathEnhancers=${pathRoadmap}/enhancer_regions/selected_enhancers.bed --pathOutput=${pathResults}/selected_promoters_enhancers_in_contact_simulated_data.txt
fi


if [ -e ${pathResults}/FOCS_promoters_enhancers_in_contact_simulated_data.txt ]; then
    echo "already done"
else
    perl ${pathScripts}/extract.regions.in.contact.pl --pathContacts=${pathInteractions}/simulations_merged.txt  --pathPromoters=${pathRoadmap}/promoter_regions/promoters_FOCS.bed  --pathEnhancers=${pathRoadmap}/enhancer_regions/enhancers_FOCS.bed --pathOutput=${pathResults}/FOCS_promoters_enhancers_in_contact_simulated_data.txt
fi

###################################################################################

