#!/bin/bash

################ Input parameters ################################

export species=$1         # i.e : human or mouse
export dataset=$2        # i.e : FANTOM5 ENCODE RoadmapEpigenomics FOCS_GRO_seq restriction_fragments

################ Export paths #####################################

export path=/beegfs/data/necsulea/RegulatoryLandscapes
export pathFinalData=/beegfs/data/necsulea/RegulatoryLandscapesManuscript/

if [ "${dataset}" = "restriction_fragments" ]; then
    export pathElements=${pathFinalData}/SupplementaryDataset1/${species}
else
    export pathElements=${pathFinalData}/SupplementaryDataset4/${species}/${dataset}
fi

export pathExons=/beegfs/data/alaverre/Regulatory_landscape/data/exons/${species}_exons_Ensembl94.txt
export pathRepeatMasker=${path}/data/RepeatMasker/${species}
export pathResults=${path}/results/sequence_composition/${species}/${dataset}
export pathScripts=${path}/scripts/sequence_composition

####################################################################

if [ "${species}" = "mouse" ]; then
    export genome=mm10
elif [ "${species}" = "human" ]; then
    export genome=hg38
fi

if [ "${dataset}" = "restriction_fragments" ]; then
    export CoordSuffix=frag_coords_${genome}.bed
    export coordConvention=1_closed_end
else
    export CoordSuffix=enhancer_coordinates.bed
    export coordConvention=0_open_end
fi

################ Create output directory #####################

if [ -e "${pathResults}" ]; then
    echo "dir output already there"
else
    mkdir -p ${pathResults}
fi

################ Running overlap with repeats #####################

echo "#!/bin/bash" > ${pathScripts}/sub_overlap_repeats

echo "perl ${pathScripts}/overlap.repeats.pl --pathCoords=${pathElements}/${CoordSuffix} --coordConvention=${coordConvention} --pathExonBlocks=${pathExons} --pathRepeatMasker=${pathRepeatMasker}/RepeatMasker_UCSC.txt.gz --pathOutput=${pathResults}/overlap_exons_repeats.txt" >> ${pathScripts}/sub_overlap_repeats

sbatch -p normal --time=1:00:00 --mem=80GB -c 1 -o ${pathScripts}/std_output_${species}_${dataset} -e ${pathScripts}/std_error_${species}_${dataset} ${pathScripts}/sub_overlap_repeats

#####################################################################
