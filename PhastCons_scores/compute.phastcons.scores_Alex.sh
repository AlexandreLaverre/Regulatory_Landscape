#!/bin/bash

################ Input parameters #####################
export species=$1
export phast=$2
export enhancer=$3
export masked_exons=$4

################ Export paths #####################
export path=/beegfs/data/alaverre/Regulatory_landscape/

export pathEnhancers=${path}/result/Supplementary_dataset3_annotations/${species}/coord_enh
export pathPhastCons=${path}/data/phastcons/${species}

export pathResults=${path}/result/Supplementary_dataset6_regulatory_landscape_evolution/${species}/enhancers_conservation/PhastCons
export pathScripts=${path}/scripts/sequence_evolution/phastcons

################ Define aligments input files #####################

if [ ${species} = "mouse" ]; then
  export aligned_species=60way
  export chromosomes=({1..19} X Y)
elif [ ${species} = "human" ]; then
  export chromosomes=({1..22} X Y)
  export aligned_species=100way
fi

if [ ${phast} = "placental" ]; then
	export suffixPhast=.phastCons${aligned_species}.placental.wigFix.gz
elif [ ${phast} = "vertebrates" ]; then
	export suffixPhast=.phastCons${aligned_species}.wigFix.gz
fi

if [ ${masked_exons} = "TRUE" ]; then
  release=94
  export suffixExons=ExonBlocks_Ensembl${release}
  export pathExons=${path}/data/ensembl_annotations/${species}/${suffixExons}.txt
else
  export suffixExons=Unmasked_Exons
  export pathExons=NA
fi

################ Create dir output #####################

if [ -e ${pathResults} ]; then
    echo "dir output already there"
else
    mkdir ${pathResults}
fi

################ Running Compute Phastcons scores #####################

for chr in ${chromosomes[@]}
do
    export pathPhast=${pathPhastCons}/chr${chr}${suffixPhast}

    if [ -e ${pathPhast} ]; then
      if [ -e ${pathResults}/PhastCons_Enhancers_chr${chr}_${enhancer}.txt ]; then
	    echo "already done"

	  else
    	echo "#!/bin/bash" > bsub_script
	    echo "#PBS -o std_output_phast.txt" >> bsub_script
	    echo "#PBS -e std_error_phast.txt" >> bsub_script
	    
	    echo "perl ${pathScripts}/compute.phastcons.scores.pl --pathCoords=${pathEnhancers}/PromoterCoords_${enhancer}.bed  --pathMaskExonBlocks=${pathExons} --pathPhastCons=${pathPhast} --chr=chr${chr} --pathOutput=${pathResults}/PhastCons_Enhancers_MaskedExons_chr${chr}_${suffixEnh}.txt" >> bsub_script

		  qsub -P P_biometr -q mc_highmem_long -l s_rss=30G,sps=1 -pe multicores 1 -o ${pathScripts}/std_output_phastcons_promoters_${species}.txt -e ${pathScripts}/std_error_phastcons_promoters_${species}.txt ${pathScripts}/bsub_script

	    fi
    fi
done

#####################################################################
