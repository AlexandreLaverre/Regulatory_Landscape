#!/bin/bash

sp1=$1 # ex: human
sp2=$2 # ex : mouse
file=$3 # ex : CAGE, merged_interacted_fragments
pathResult="/beegfs/data/alaverre/Regulatory_landscape/result/genome_alignment/${sp1}2other/${sp1}2${sp2}/"
pathScript="/beegfs/data/alaverre/Regulatory_landscape/script/align_sequences/"

for list in `ls ${pathResult}list/${file}* | xargs -n 1 basename` 
do
echo "#!/bin/bash" > ${pathResult}/pecan_alignments/pecan_${sp1}2${sp2}_${list} 
echo "#SBATCH -o ${pathResult}/pecan_alignments/std_output_${list}.txt" >> ${pathResult}/pecan_alignments/pecan_${sp1}2${sp2}_${list}
echo "#SBATCH -e ${pathResult}/pecan_alignments/std_error_${list}.txt" >> ${pathResult}/pecan_alignments/pecan_${sp1}2${sp2}_${list}	
echo "source /beegfs/home/alaverre/.bashrc" >> ${pathResult}/pecan_alignments/pecan_${sp1}2${sp2}_${list}
echo "python ${pathScript}run_pecan.py ${list} ${sp1} ${sp2} ${file}" >> ${pathResult}/pecan_alignments/pecan_${sp1}2${sp2}_${list}
sbatch -p normal --time=24:00:00 --mem=10GB -c 1 ${pathResult}/pecan_alignments/pecan_${sp1}2${sp2}_${list}
done

