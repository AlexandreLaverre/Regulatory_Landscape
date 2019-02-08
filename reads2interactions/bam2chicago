#!/bin/bash
specie=$1
genome=$2
author_ctype=$3

path_bam="/beegfs/data/alaverre/result/HiCup/${specie}/${author_ctype}/"
path_data="/beegfs/data/alaverre/data/CHICAGO_files/${specie}/${genome}/"
path_result="/beegfs/data/alaverre/result/CHICAGO/${specie}/${author_ctype}/"
mkdir -p ${path_result}
for sample in `ls ${path_bam}*.dedup.bam`
do
echant=$(echo ${sample} | rev | cut -d / -f 1  | rev | cut -d . -f 1)
echo "#!/bin/bash" > ${path_result}sub_bam2chic_${echant}
echo "#SBATCH -o ${path_result}std_output_bam2chic_${echant}.txt" >> ${path_result}sub_bam2chic_${echant}
echo "#SBATCH -e ${path_result}std_error_bam2chic_${echant}.txt" >> ${path_result}sub_bam2chic_${echant}	
echo "source /beegfs/home/alaverre/.bashrc" >> ${path_result}sub_bam2chic_${echant} 
echo "bam2chicago.sh ${sample} ${path_data}Digest_${genome}_HindIII_None.txt.baitmap ${path_data}Digest_${genome}_HindIII_None.txt.rmap ${path_result}${echant}" >> ${path_result}sub_bam2chic_${echant}
sbatch -p normal --time 72:00:00 --mem=1G -c 1 ${path_result}sub_bam2chic_${echant}
done

