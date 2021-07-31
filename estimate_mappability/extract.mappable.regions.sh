#/bin/bash

export sp=$1
export readlen=$2
export cluster=$3

##############################################################

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/${USER}/RegulatoryLandscapes
fi

export pathGenomeIndexes=${path}/data/genome_sequences/${sp}
export pathResults=${path}/results/mappability/readlength${readlen}/${sp}
export pathScripts=${path}/scripts/estimate_mappability

##############################################################

for chr in {1..32} X Y Z W 
do
    
    if [ -e ${pathResults}/fake_reads_chr${chr}_mappedregions.txt ]; then
	echo "already done"
    else
	if [ -e ${pathResults}/fake_reads_posinfo_chr${chr}.fa.gz ]&&[ -e ${pathResults}/fake_reads_chr${chr}.sam.gz ]; then

	    echo "#!/bin/bash" >  ${pathScripts}/bsub_script_estmap
	    
	    if [ ${cluster} = "pbil" ]; then
		echo "#SBATCH --job-name=estmap_${sp}_${chr}_${readlen}" >>  ${pathScripts}/bsub_script_estmap
		echo "#SBATCH --partition=normal" >>  ${pathScripts}/bsub_script_estmap
		echo "#SBATCH --output=${pathScripts}/std_out_estmap_${sp}_${chr}" >>  ${pathScripts}/bsub_script_estmap
		echo "#SBATCH --error=${pathScripts}/std_err_estmap_${sp}_${chr}" >>  ${pathScripts}/bsub_script_estmap
		echo "#SBATCH --cpus-per-task=1" >>  ${pathScripts}/bsub_script_estmap ## 8 CPUs
		echo "#SBATCH --time=8:00:00" >>  ${pathScripts}/bsub_script_estmap ## 8 hours
		echo "#SBATCH --mem=10G" >>  ${pathScripts}/bsub_script_estmap ## 4g per CPU
	    fi
	    
	    echo "perl ${pathScripts}/extract.mappable.regions.pl --pathReads=${pathResults}/fake_reads_posinfo_chr${chr}.fa.gz --chr=${chr} --pathAlignment=${pathResults}/fake_reads_chr${chr}.sam.gz --readLength=${readlen} --pathOutput=${pathResults}/fake_reads_chr${chr}_mappedregions.txt" >> bsub_script_estmap
	    
	    if [ ${cluster} = "in2p3" ]; then
		qsub -q huge -l s_rss=10G,sps=1 -o ${pathScripts}/std_output_${readlen}_${chr}.txt -e ${pathScripts}/std_error_${readlen}_${chr}.txt ${pathScripts}/bsub_script_estmap
	    fi
	    
	    if [ ${cluster} = "cloud" ]; then
		chmod a+x ${pathScripts}/bsub_script_estmap
		${pathScripts}/bsub_script_estmap
	    fi
	    
	    if [ ${cluster} = "pbil" ]; then
		sbatch ${pathScripts}/bsub_script_estmap
	    fi
	fi
    fi

done


##############################################################
