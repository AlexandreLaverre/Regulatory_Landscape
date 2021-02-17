#!/bin/bash

#######################################################################################

export pathOrigin=/beegfs/data/alaverre/Regulatory_landscape/result/Supplementary_dataset6_regulatory_landscape_evolution
export pathDest=/beegfs/data/necsulea/RegulatoryLandscapesManuscript/SupplementaryDataset7

#######################################################################################

for ref in human mouse
do
    for enh in ENCODE FANTOM5 FOCS_GRO_seq RoadmapEpigenomics
    do
	if [ ${enh} = "FANTOM_GRO_seq" ]; then
	    export origin_enh=GRO_seq
	fi

	if [ ${enh} = "RoadmapEpigenomics" ]; then
	    export origin_enh=Roadmap
	fi

	if [ ${enh} = "ENCODE" ]; then
	    export origin_enh=ENCODE
	fi
	
	if [ ${enh} = "FANTOM5" ]; then
	    export origin_enh=CAGE
	fi

	
	for tg in human macaque mouse rat rabbit dog cow elephant opossum chicken
	do
	    
	    if [ -e ${pathOrigin}/${ref}/synteny_conservation/${origin_enh}/${ref}2${tg}_${dataset}_synteny.txt_unique ]; then
		echo ${pathOrigin}/${ref}/synteny_conservation/${origin_enh}/${ref}2${tg}_${dataset}_synteny.txt_unique 
	    else
		echo "not found "${ref} ${tg} ${enh} {dataset}
	    fi

	done
    done
done

##############################################################################################################
