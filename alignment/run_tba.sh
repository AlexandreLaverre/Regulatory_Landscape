#!/bin/bash

pathResult="/beegfs/data/alaverre/result/genome_alignment/mouse2other/"
pathScript="/beegfs/data/alaverre/script/"

for file in `ls ${path}xaa | xargs -n 1 basename`
do
echo "#!/bin/bash" > ${pathResult}tba_${file} 
echo "#SBATCH -o ${pathResult}std_output_${file}.txt" >> ${pathResult}tba_${file}
echo "#SBATCH -e ${pathResult}std_error_${file}.txt" >> ${pathResult}tba_${file}	
echo "source /beegfs/home/alaverre/.bashrc" >> ${pathResult}tba_${file}
echo "/usr/bin/X11/python ${pathScript}run_tba.py ${file}" >> ${pathResult}tba_${file}
sbatch -p normal --time=4:00:00 --mem=10GB -c 1 ${pathResult}tba_${file}
done

