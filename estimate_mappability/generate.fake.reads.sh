#/bin/bash

export sp=$1
export readlen=$2
export cluster=$3

##############################################################

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/${USER}/RegulatoryLandscapes
fi

export pathSequence=${path}/data/genome_sequences/${sp}
export pathResults=${path}/results/mappability/readlength${readlen}/${sp}
export pathScripts=${path}/scripts/estimate_mappability

export release=94

##############################################################

if [ -e ${pathResults} ]; then 
    echo "path output already there"
else
    mkdir ${pathResults}
fi

##############################################################

if [ -e ${pathSequence}/genome_Ensembl${release}.fa ]; then
    export suffix=fa
else
    if [ -e ${pathSequence}/genome_Ensembl${release}.fa.gz ]; then
	export suffix=fa.gz
    else
	echo "cannot find fasta file"
	exit
    fi
fi

##############################################################

for chr in {1..20} X Y 
do
    
    echo "#!/bin/bash" >  ${pathScripts}/bsub_script_estmap

    if [ ${cluster} = "pbil" ]; then
	echo "#SBATCH --job-name=estmap_${sp}_${chr}" >>  ${pathScripts}/bsub_script_estmap
	echo "#SBATCH --partition=normal" >>  ${pathScripts}/bsub_script_estmap
	echo "#SBATCH --output=${pathScripts}/std_out_${sp}_${chr}_${readlen}" >>  ${pathScripts}/bsub_script_estmap
	echo "#SBATCH --error=${pathScripts}/std_err_${sp}_${chr}_${readlen}" >>  ${pathScripts}/bsub_script_estmap
	echo "#SBATCH --cpus-per-task=1" >>  ${pathScripts}/bsub_script_estmap ## 1 CPU
	echo "#SBATCH --time=5:00:00" >>  ${pathScripts}/bsub_script_estmap ## 5 hours
	echo "#SBATCH --mem=5G" >>  ${pathScripts}/bsub_script_estmap ## 4g per CPU
    fi
    
    echo "perl ${pathScripts}/generate.fake.reads.pl --pathGenome=${pathSequence}/genome_Ensembl${release}.${suffix} --readLength=${readlen} --step=5 --chr=${chr} --pathOutput=${pathResults}/fake_reads_posinfo_chr${chr}.fa" >> ${pathScripts}/bsub_script_estmap

    echo "gzip ${pathResults}/fake_reads_posinfo_chr${chr}.fa" >> ${pathScripts}/bsub_script_estmap
    
    if [ ${cluster} = "pbil" ]; then
	sbatch ${pathScripts}/bsub_script_estmap
    fi

done

##############################################################



