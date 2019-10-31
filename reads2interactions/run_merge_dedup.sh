#!/bin/bash
specie=$1			#ex : mouse
author_ctype=$2			#ex : Schoenfelder/ESC
name=$3				#ex : ESC_rep1

path_result="/beegfs/data/alaverre/result/HiCup/${specie}/${author_ctype}/"

files=`ls ${path_result}*.hicup.bam`

echo "#!/bin/bash" > ${path_result}sub_dedup_${name}
echo "#SBATCH -o ${path_result}std_output_${name}.txt" >> ${path_result}sub_dedup_${name}
echo "#SBATCH -e ${path_result}std_error_${name}.txt" >> ${path_result}sub_dedup_${name}
echo "source /beegfs/home/alaverre/.bashrc" >> ${path_result}sub_dedup_${name}

echo "samtools merge -n ${path_result}${name}_sorted.bam `echo ${files}` -@ 10" >> ${path_result}sub_dedup_${name}
echo "hicup_deduplicator --zip ${path_result}${name}_sorted.bam  --outdir ${path_result} --threads 10" >> ${path_result}sub_dedup_${name}

sbatch -p normal --time=24:00:00 --mem=10G -c 10 --nodelist pbil-deb[23-29] ${path_result}sub_dedup_${name}
