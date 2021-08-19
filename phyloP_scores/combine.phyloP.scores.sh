#!/bin/bash

################ Input parameters #####################
export species=$1         # i.e : human or mouse
export score=$2		  # i.e : phyloP or phastCons
export way=$3             # i.e : 30 60 or 100way
export enhancer=$4        # i.e : FANTOM5 ENCODE RoadmapEpigenomics FOCS_GRO_seq restriction_fragments
export masked_exons=$5    # i.e : TRUE or FALSE

################ Export paths #####################
export path=/beegfs/data/alaverre/Regulatory_landscape/
export pathAnouk=/beegfs/data/necsulea/RegulatoryLandscapesManuscript/


if [ "${enhancer}" = "restriction_fragments" ]; then
export pathEnhancers=${pathAnouk}/SupplementaryDataset1/${species}
else
export pathEnhancers=${pathAnouk}/SupplementaryDataset4/${species}/${enhancer}
fi

export pathphyloP=${path}/data/${score}/${species}/${way}

if [ ${USER} = "alaverre" ]; then
    export pathResults=${path}/result/${score}/${species}/${enhancer}
    export pathScripts=${path}/scripts/phyloP_scores
fi

if [ ${USER} = "necsulea" ]; then
    export pathResults=/beegfs/data/necsulea/RegulatoryLandscapes/results/${score}/${species}/${enhancer}
    export pathScripts=/beegfs/data/necsulea/RegulatoryLandscapes/scripts/phyloP_scores
fi

################################################################################

if [ "${species}" = "mouse" ]; then
  export chromosomes=({1..19} X Y)
elif [ "${species}" = "human" ]; then
  export chromosomes=({1..22} X Y)
fi

################################################################################

if [ "${masked_exons}" = "TRUE" ]; then
    export suffixExons=MaskedExons_Ensembl94
else
    export suffixExons=Unmasked
fi

################################################################################

if [ -e ${pathResults}/${score}_${way}_${suffixExons}.txt ]; then
    echo "path output already exists, not doing anything"
    exit
fi

################################################################################

for chr in "${chromosomes[@]}"
do
    if [ -e ${pathResults}/${score}_${way}_${suffixExons}.txt ]; then
	sed '1d' ${pathResults}/${score}_${way}_chr${chr}_${suffixExons}.txt >> ${pathResults}/${score}_${way}_${suffixExons}.txt
    else
	cat ${pathResults}/${score}_${way}_chr${chr}_${suffixExons}.txt > ${pathResults}/${score}_${way}_${suffixExons}.txt
    fi
done

################################################################################
