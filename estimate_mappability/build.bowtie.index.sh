#/bin/bash

export sp=$1
export cluster=$2

##############################################################

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/${USER}/RegulatoryLandscapes
fi

export pathGenomeIndexes=${path}/data/genome_sequences/${sp}
export pathScripts=${path}/scripts/estimate_mappability

export release=94

## bowtie2 version: 2.3.4.3

##############################################################

echo "#!/bin/bash" >  ${pathScripts}/bsub_script_bowtie_index

if [ ${cluster} = "pbil" ]; then
    echo "#SBATCH --job-name=bowtie_index_${sp}_${chr}" >>  ${pathScripts}/bsub_script_bowtie_index
    echo "#SBATCH --partition=normal" >>  ${pathScripts}/bsub_script_bowtie_index
    echo "#SBATCH --output=${pathScripts}/std_out_index_${sp}" >>  ${pathScripts}/bsub_script_bowtie_index
    echo "#SBATCH --error=${pathScripts}/std_err_index_${sp}" >>  ${pathScripts}/bsub_script_bowtie_index
    echo "#SBATCH --cpus-per-task=8" >>  ${pathScripts}/bsub_script_bowtie_index ## 8 CPUs
    echo "#SBATCH --time=8:00:00" >>  ${pathScripts}/bsub_script_bowtie_index ## 8 hours
    echo "#SBATCH --mem=4G" >>  ${pathScripts}/bsub_script_bowtie_index ## 4g per CPU
fi

echo "bowtie2-build --threads 8 ${pathGenomeIndexes}/genome_Ensembl${release}.fa ${pathGenomeIndexes}/genome_Ensembl${release}" >>  ${pathScripts}/bsub_script_bowtie_index

if [ ${cluster} = "pbil" ]; then
    sbatch ${pathScripts}/bsub_script_bowtie_index
fi

##############################################################
