#!/bin/bash

sp_origin=$1
sp_target=$2
data=$3
nb_part=$4
pathScript=/beegfs/data/alaverre/Regulatory_landscape/test/script

source ~/.bashrc
conda activate /beegfs/data/alaverre/Tools/envs
snakemake --unlock 
snakemake -s ${pathScript}/Snakefile -j 500 --config sp_origin=${sp_origin} sp_target=${sp_target} data=${data} nb_part=${nb_part} --rerun-incomplete --cluster "sbatch --job-name=snakejob.{rule}.${sp_origin}2${sp_target}.${data} -p normal -N 1 -o /beegfs/data/alaverre/Regulatory_landscape/test/result/${sp_origin}2other/${sp_origin}2${sp_target}/log_files/slurm.out_${data} -e /beegfs/data/alaverre/Regulatory_landscape/test/result/${sp_origin}2other/${sp_origin}2${sp_target}/log_files/slurm.err_${data} -c {params.threads} --mem={params.mem} -t {params.time}"
echo "Done!"
sleep 1h
