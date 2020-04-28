#!/bin/bash

################ Input parameters #####################
export species=$1         # i.e : human or mouse
export phast=$2           # i.e : placental or vertebrates
export enhancer=$3        # i.e : CAGE ENCODE RoadMap GRO_seq
export masked_exons=$4    # i.e : TRUE or FALSE

################ Export paths #####################
export path=/beegfs/data/alaverre/Regulatory_landscape

export pathEnhancers=${path}/result/Supplementary_dataset3_annotations/${species}/coord_enh
export pathPhastCons=${path}/data/phastcons/${species}

export pathResults=${path}/result/Supplementary_dataset6_regulatory_landscape_evolution/${species}/enhancers_conservation/PhastCons
export pathScripts=${path}/script/PhastCons_scores

################ Define aligments input files #####################

if [ "${species}" = "mouse" ]; then
  export aligned_species=60way
  export chromosomes=({1..19} X Y)
  export genome=mm10
elif [ "${species}" = "human" ]; then
  export chromosomes=({1..22} X Y)
  export aligned_species=30way
  export genome=hg38
fi

if [ "${phast}" = "placental" ]; then
	export suffixPhast=.phastCons${aligned_species}.placental.wigFix.gz
elif [ "${phast}" = "vertebrates" ]; then
	export suffixPhast=.phastCons${aligned_species}.wigFix.gz
fi

if [ "${masked_exons}" = "TRUE" ]; then
  release=94
  export suffixExons=MaskedExons_Ensembl${release}
  export pathExons=${path}/test/data/exons/${species}_all_exons_merged.bed
else
  export suffixExons=Unmasked_Exons
  export pathExons=NA
fi

################ Create dir output #####################

if [ -e "${pathResults}" ]; then
    echo "dir output already there"
else
    mkdir "${pathResults}"
fi

################ Running Compute Phastcons scores #####################

for chr in "${chromosomes[@]}"
do
    export pathPhast=${pathPhastCons}/chr${chr}${suffixPhast}

    if [ -e "${pathPhast}" ]; then
      if [ -e "${pathResults}"/PhastCons_Enhancers_chr"${chr}"_"${enhancer}".txt ]; then
	    echo "already done"

	  else
    	echo "#!/bin/bash" > PhastCons_${phast}_${enhancer}_chr${chr}_${suffixExons}
	    echo "#PBS -o std_output_phast.txt" >> PhastCons_${phast}_${enhancer}_chr${chr}_${suffixExons}
	    echo "#PBS -e std_error_phast.txt" >> PhastCons_${phast}_${enhancer}_chr${chr}_${suffixExons}
	    
	    echo "perl ${pathScripts}/compute.phastcons.scores.pl --pathCoords=${pathEnhancers}/${enhancer}_enhancer_genomic_positions_${genome}.bed  --pathMaskExonBlocks=${pathExons} --pathPhastCons=${pathPhast} --chr=chr${chr} --pathOutput=${pathResults}/PhastCons_${phast}_${aligned_species}_${enhancer}_chr${chr}_${suffixExons}.txt" >> PhastCons_${phast}_${enhancer}_chr${chr}_${suffixExons}

      bash ${pathScripts}/PhastCons_${phast}_${enhancer}_chr${chr}_${suffixExons}
      #sbatch -p normal --exclude=pbil-deb11 --time=24:00:00 --mem=30GB -c 1 -o ${pathScripts}/std_output_phastcons_enhancers_${species}.txt -e ${pathScripts}/std_error_phastcons_enhancers_${species}.txt ${pathScripts}/PhastCons_${phast}_${enhancer}_chr${chr}_${suffixExons}
	    fi
    fi
done

#####################################################################
