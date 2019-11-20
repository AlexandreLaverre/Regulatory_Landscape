#!/bin/bash

export sample=$1
export type=$2
export cluster=$3

###################################################################################

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/necsulea/RegulatoryLandscapes
else
    if [ ${cluster} = "in2p3" ]; then
	export path=/sps/biometr/necsulea/RegulatoryLandscapes
    else
	echo "unknown cluster"
	exit
    fi
fi

export pathRoadmap=${path}/data/RoadmapEpigenomics/hg19
export pathCoverage=${pathRoadmap}/signal_tracks
export pathEnhancers=${pathRoadmap}/regulatory_regions/enhancer_regions
export pathPromoters=${pathRoadmap}/regulatory_regions/promoter_regions
export pathScripts=${path}/scripts/process_RoadmapEpigenomics_enhancers

if [ ${type} = "enhancers" ]; then
    export pathResults=${pathEnhancers}
    export prefix="all_enhancers"
else
    if [ ${type} = "promoters" ]; then
	export pathResults=${pathPromoters}
	export prefix="all_promoters"
    else
	echo "unknown type"
	exit
    fi
fi

###################################################################################

if [ -e ${pathResults}/${prefix}_coverage_${sample}.txt ]; then
    echo "already done"
else
    if [ -e ${pathCoverage}/${sample}.bedGraph.gz ]; then
	
	echo "#!/bin/bash" >  ${pathScripts}/bsub_script_coverage
	
	echo "perl ${pathScripts}/compute.coverage.pl --pathCoordinates=${pathResults}/${prefix}.bed  --pathCoverage=${pathCoverage}/${sample}.bedGraph.gz --pathOutput=${pathResults}/${prefix}_coverage_${sample}.txt" >> ${pathScripts}/bsub_script_coverage
    	
	
	if [ ${cluster} = "in2p3" ]; then
	    qsub -q huge -l s_rss=6G,sps=1 -o ${pathScripts}/std_output_coverage_${sample}_${type}.txt -e ${pathScripts}/std_error_coverage_${sample}_${type}.txt ${pathScripts}/bsub_script_coverage
	fi
    fi	    
done

###################################################################################
