#!/bin/bash
specie=$1
author=$2

path_data="/pandata/alaverre/data/reads/sra"
path_result="/pandata/alaverre/data/reads_PC_HIC/${specie}/${author}"

for file in `ls ${path_data}/*.sra | xargs -n 1 basename`
do
mkdir -p ${path_result}/${file}
echo "#!/bin/bash" > ${path_result}/bsub_fastq_${file}
echo "#PBS -o ${path_result}/std_output.txt" >> ${path_result}/bsub_fastq_${file}
echo "#PBS -e ${path_result}/std_error.txt" >> ${path_result}/bsub_fastq_${file}
echo "source /panhome/alaverre/.bashrc" >> ${path_result}/bsub_fastq_${file}
echo "parallel-fastq-dump --split-files --gzip --threads 10 --sra-id ${path_data}/${file} --outdir ${path_result}/${file} --tmpdir ${path_result}/${file}" >> ${path_result}/bsub_fastq_${file}
qsub -q q1day -l nodes=1:ppn=10 ${path_result}/bsub_fastq_${file}
done

