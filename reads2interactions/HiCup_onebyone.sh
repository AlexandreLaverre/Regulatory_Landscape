#!/bin/bash
specie=$1
genome=$2
author_ctype=$3
path_data="/beegfs/data/alaverre/data/reads_PC_HIC/${specie}/${author_ctype}/"
path_result="/beegfs/data/alaverre/result/HiCup/${specie}/${author_ctype}/"

for sample in `ls ${path_data}*1.fastq.gz | xargs -n 1 basename | cut -d : -f 1`
do
mkdir -p ${path_result}${sample}
## Config files for HiCUP ##
echo "Outdir: ${path_result}${sample}/" > ${path_result}${sample}/config_files.txt
echo "Zip: 1" >> ${path_result}${sample}/config_files.txt
echo "Keep: 0" >> ${path_result}${sample}/config_files.txt
echo "Threads: 10" >> ${path_result}${sample}/config_files.txt
echo "Bowtie2: /beegfs/home/alaverre/Tools/bowtie2/bowtie2" >> ${path_result}${sample}/config_files.txt 
echo "Digest: /beegfs/data/alaverre/data/genome_ref/Digest_${genome}_HindIII_None.txt" >> ${path_result}${sample}/config_files.txt
echo "Index: /beegfs/data/alaverre/data/genome_ref/${genome}" >> ${path_result}${sample}/config_files.txt
echo "R: /beegfs/data/soft/R-3.5.1/bin/R" >> ${path_result}${sample}/config_files.txt
echo "Shortest: 150" >> ${path_result}${sample}/config_files.txt
echo "Longest: 800" >> ${path_result}${sample}/config_files.txt

echo "${path_data}${sample}:1.fastq.gz" >> ${path_result}${sample}/config_files.txt
echo "${path_data}${sample}:2.fastq.gz" >> ${path_result}${sample}/config_files.txt
echo "" >> ${path_result}${sample}/config_files.txt
sed -i '$d' ${path_result}${sample}/config_files.txt

## HiCUP ##
echo "#!/bin/bash" > ${path_result}${sample}/bsub_hicup_${sample} 
echo "#SBATCH -o ${path_result}${sample}/std_output.txt" >> ${path_result}${sample}/bsub_hicup_${sample}
echo "#SBATCH -e ${path_result}${sample}/std_error.txt" >> ${path_result}${sample}/bsub_hicup_${sample}	
echo "source /beegfs/home/alaverre/.bashrc" >> ${path_result}${sample}/bsub_hicup_${sample}
echo "hicup --config  ${path_result}${sample}/config_files.txt" >> ${path_result}${sample}/bsub_hicup_${sample}
sbatch -p normal --time=30:00:00 --mem=8GB -c 10 -C "haswell" ${path_result}${sample}/bsub_hicup_${sample}
done
