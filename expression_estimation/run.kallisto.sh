#!/bin/bash

export species=$1
export type=$2
export cluster=$3

####################################################################################

if [ ${cluster} = "pbil" ]; then
    export path=/bgeefs/data/necsulea/RegulatoryLandscapes
fi

if [ ${cluster} = "in2p3" ]; then
    export path=/sps/biometr/necsulea/RegulatoryLandscapes
fi

if [ ${cluster} = "cloud" ]; then
    export path=/mnt/RegulatoryLandscapes
fi

####################################################################################

export pathDocs=${path}/docs
export pathRNASeq=${path}/data/RNASeq/${species}
export pathResults=${path}/results/expression_estimation/${species}
export pathIndexes=${path}/results/kallisto_indexes/${species}
export pathScripts=${path}/scripts/expression_estimation

export release=94

####################################################################################

for sample in `cut -f 1 ${pathDocs}/RNASeq_Samples_${species}.txt | uniq | grep -v SampleID`
do
    export pathRead1=""
    export pathRead2=""

    export library="single"
    
    for file in `ls ${pathRNASeq} | grep ^${sample} | grep _2.fastq.gz`
    do
	export library="paired"
    done
        	
    echo ${sample} ${library}
    
    if [ -e ${pathResults}/${sample}/kallisto_${type}/abundance.tsv ]; then
	echo "already done"
    else

	if [ -e ${pathResults}/${sample}/kallisto_${type} ]; then
	    echo "output dir already there"
	else
	    mkdir -p ${pathResults}/${sample}/kallisto_${type}
	fi

	

	if [ ${cluster} = "cloud" ]; then
	    if [ ${library} = "single" ]; then
		kallisto quant --single -l 200.0 -s 20 --bias -t 4 -o ${pathResults}/${sample}/kallisto_${type} --index ${pathIndexes}/${type} <(zcat ${pathRNASeq}/${sample}_*_1.fastq.gz)
	    fi

	    if [ ${library} = "paired" ]; then
		kallisto quant -l 200.0 -s 20 --bias -t 4 -o ${pathResults}/${sample}/kallisto_${type} --index ${pathIndexes}/${type} <(zcat ${pathRNASeq}/${sample}_*_1.fastq.gz) <(zcat ${pathRNASeq}/${sample}_*_2.fastq.gz)
	    fi

	fi
    fi
done

####################################################################################
