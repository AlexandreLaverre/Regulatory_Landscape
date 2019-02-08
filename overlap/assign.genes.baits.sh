#!/bin/bash



#####################################################################



export path=/pandata/alaverre

export pathPromoterCaptureHiC=${path}/data
export pathAnnotations=${path}/data/gene_annotations/human

export pathScripts=${path}/scripts/annotate_chromatin_contacts



export release=91



#####################################################################



# perl ${pathScripts}/assign.genes.baits.pl --pathPromoterCaptureHiC="to_be_defined" --pathTranscriptInfo=${pathAnnotations}/mouse/TranscriptInfo_Ensembl${release}.txt --pathGeneNames=${pathAnnotations}/mouse/GeneNames_Ensembl${release}.txt --maxDistance=5000 --pathOutput="to_be_defined"



perl ./assign.genes.baits.pl --pathPromoterCaptureHiC=${pathPromoterCaptureHiC}/human_promoter-other_hg38_annot5kb_filtred.txt  --pathTranscriptInfo=${pathAnnotations}/TranscriptInfo_Filtered_Ensembl${release}.txt --pathGeneNames=${pathAnnotations}/GeneNames_Ensembl${release}.txt --maxDistance=5000 --pathOutput=${pathPromoterCaptureHiC}/human_promoter-other_hg38_annot5kb_filtred_annotreal.txt


###################################################################

