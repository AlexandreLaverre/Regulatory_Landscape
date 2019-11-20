#!/bin/bash

#####################################################################

export path=/sps/biometr/necsulea/RegulatoryLandscapes
export pathEnhancers=${path}/data/FANTOM5
export pathAlignments=${path}/data/liftOver_files
export pathScripts=${path}/scripts/process_FANTOM5_enhancers

#####################################################################

for genome in mm9 hg19
do
    if [ ${genome} = "mm9" ]; then
	export sp="mouse"
	export assembly="mm10"
	export liftOverFile=mm9ToMm10.over.chain.gz
    fi

    if [ ${genome} = "hg19" ]; then
	export sp="human"
	export assembly="hg38"
	export liftOverFile=hg19ToHg38.over.chain.gz
    fi
    
    liftOver ${pathEnhancers}/${genome}/${sp}_permissive_enhancers_phase_1_and_2.bed ${pathAlignments}/${liftOverFile} ${pathEnhancers}/${genome}/${sp}_lifted_to_${assembly}_permissive_enhancers_phase_1_and_2.bed ${pathEnhancers}/${genome}/${sp}_lifted_to_${assembly}_permissive_enhancers_phase_1_and_2.unmapped
   
done

#####################################################################
