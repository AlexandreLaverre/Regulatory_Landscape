#!/bin/bash

path="/beegfs/data/alaverre/result/genome_alignment/mouse2other/"

for file in `ls ${path}xaa`
do
echo "#!/bin/bash" > ${path}tba_${file} 
echo "#SBATCH -o ${path}std_output_${file}.txt" >> ${path}tba_${file}
echo "#SBATCH -e ${path}std_error_${file}.txt" >> ${path}tba_${file}	
echo "source /beegfs/home/alaverre/.bashrc" >> ${path}tba_${file}
echo "python ${path}/run_tba.py ${file}" >> ${path}tba_${file}
sbatch -p normal --time=4:00:00 --mem=10GB -c 1 ${path}tba_${file}
