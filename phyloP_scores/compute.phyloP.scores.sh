#!/bin/bash

################ Input parameters #####################
export species=$1         # i.e : human or mouse
export way=$2             # i.e : 30 60 or 100way
export enhancer=$3        # i.e : FANTOM5 ENCODE RoadmapEpigenomics FOCS_GRO_seq restriction_fragments

################ Export paths #####################
export path=/beegfs/data/alaverre/Regulatory_landscape/
export pathAnouk=/beegfs/data/necsulea/RegulatoryLandscapesManuscript/


if [ "${enhancer}" = "restriction_fragments" ]; then
export pathEnhancers=${pathAnouk}/SupplementaryDataset1/${species}
else
export pathEnhancers=${pathAnouk}/SupplementaryDataset4/${species}/${enhancer}
fi

export pathphyloP=${path}/data/phyloP/${species}/${way}
export pathResults=${path}/result/phyloP/${species}/${enhancer}
export pathScripts=${path}/scripts/phyloP_scores

################ Define aligments input files #####################

export suffixPhylo=phyloP${way}.wigFix.gz
export suffixExons=MaskedExons_Ensembl94
export pathExons=${path}/Snakemake_folder/data/exons/${species}_all_exons_merged.bed

if [ "${species}" = "mouse" ]; then
  export chromosomes=({1..19} X Y)
  export genome=mm10
elif [ "${species}" = "human" ]; then
  export chromosomes=({1..22} X Y)
  export genome=hg38
fi


if [ "${enhancer}" = "restriction_fragments" ]; then
export CoordSuffix=frag_coords_${genome}.bed
else
export CoordSuffix=enhancer_coordinates.bed
fi

################ Create dir output #####################

if [ -e "${pathResults}" ]; then
    echo "dir output already there"
else
    mkdir "${pathResults}"
fi

################ Running Compute phyloP scores #####################

for chr in "${chromosomes[@]}"
do
    export pathPhyloP_scores=${pathphyloP}/chr${chr}.${suffixPhylo}

    if [ -e "${pathPhyloP_scores}" ]; then
      if [ -e "${pathResults}"/phyloP_${way}_chr${chr}_${suffixExons}.txt ]; then
	    echo "already done"

      else
    	    echo "#!/bin/bash" > ${pathScripts}/log/phyloP_${way}_${species}_${enhancer}_chr${chr}_${suffixExons}
	    echo "#PBS -o std_output_phyloP.txt" >> ${pathScripts}/log/phyloP_${way}_${species}_${enhancer}_chr${chr}_${suffixExons}
	    echo "#PBS -e std_error_phyloP.txt" >> ${pathScripts}/log/phyloP_${way}_${species}_${enhancer}_chr${chr}_${suffixExons}
	    
	    echo "perl ${pathScripts}/compute.phyloP.scores.pl --pathCoords=${pathEnhancers}/${CoordSuffix}   --pathMaskExonBlocks=${pathExons} --pathPhastCons=${pathPhyloP_scores} --chr=chr${chr} --pathOutput=${pathResults}/phyloP_${way}_chr${chr}_${suffixExons}.txt" >> ${pathScripts}/log/phyloP_${way}_${species}_${enhancer}_chr${chr}_${suffixExons}

      #bash ${pathScripts}/log/phyloP_${way}_${species}_${enhancer}_chr${chr}_${suffixExons}
      sbatch -p normal --time=1:00:00 --mem=40GB -c 1 -o ${pathScripts}/log/std_output_phyloP_${way}_${species}_${enhancer}_chr${chr}_${suffixExons} -e ${pathScripts}/log/std_error_phyloP_${way}_${species}_${enhancer}_chr${chr}_${suffixExons} ${pathScripts}/log/phyloP_${way}_${species}_${enhancer}_chr${chr}_${suffixExons}
       fi
    fi
done

#####################################################################
