#/bin/bash

export sp=$1
export readlen=$2
export cluster=$3

##############################################################

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/${USER}/RegulatoryLandscapes
fi

export pathGenomeIndexes=${path}/data/genome_sequences/${ref}
export pathResults=${path}/results/mappability/readlength${readlen}/${sp}
export pathScripts=${path}/scripts/estimate_mappability

export release=94

##############################################################

export pathIndex=${pathGenomeIndexes}/genome_Ensembl${release}

##############################################################

for chr in {1..32} X Y Z W
do
    
    if [ -e ${pathResults}/fake_reads_chr${chr}.sam.gz ]; then
	echo "already done"
    else
	if [ -e ${pathResults}/fake_reads_posinfo_chr${chr}.fa.gz ]; then
	    echo "#!/bin/bash" >  ${pathScripts}/bsub_script_map

	     if [ ${cluster} = "pbil" ]; then
		 echo "#SBATCH --partition=normal" >>  ${pathScripts}/bsub_script_map
		 echo "#SBATCH --output=${pathScripts}/std_out_map_${sp}_${chr}" >>  ${pathScripts}/bsub_script_map
		 echo "#SBATCH --error=${pathScripts}/std_err_map_${sp}_${chr}" >>  ${pathScripts}/bsub_script_map
		 echo "#SBATCH --cpus-per-task=8" >>  ${pathScripts}/bsub_script_map ## 8 CPUs
		 echo "#SBATCH --time=8:00:00" >>  ${pathScripts}/bsub_script_map ## 8 hours
		 echo "#SBATCH --mem=4G" >>  ${pathScripts}/bsub_script_map ## 4g per CPU
	     fi
	     
	     echo "bowtie2-align-s --wrapper basic-0 --very-sensitive -x ${pathIndex} --threads 8 --reorder --passthrough -S ${pathResults}/fake_reads_chr${chr}.sam -U ${pathResults}/fake_reads_posinfo_chr${chr}.fa.gz -f" >> ${pathScripts}/bsub_script_map
	     
	    echo "gzip ${pathResults}/fake_reads_chr${chr}.sam" >> ${pathScripts}/bsub_script_map

	    if [ ${cluster} = "pbil" ]; then
		sbatch ${pathScripts}/bsub_script_bowtie_index
	    fi
	    
	fi
    fi
done

##############################################################


