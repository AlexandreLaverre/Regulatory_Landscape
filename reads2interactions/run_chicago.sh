#!/bin/bash
sp=$1			# ex : human, mouse
genome=$2		# ex : hg38, mm10
author_ctype=$3		# ex : Schoenfelder_ESC, Mifsud_TFC
name=$4			# ex : ESC_rep1

path_data="/beegfs/data/alaverre/data/CHICAGO_files/${sp}/${genome}/"
path_result="/beegfs/data/alaverre/result/CHICAGO/${sp}/${author_ctype}/"

sample=$(ls ${path_result}*.chinput | sed 's/\t/;/g')
files=$(echo ${sample} | sed 's/ /,/g')
run=$(echo ${author_ctype} | cut -d / -f 2)

echo "#!/bin/bash" > ${path_result}chicago_${run}
echo "#SBATCH -o ${path_result}std_output_chicago.txt" >> ${path_result}chicago_${run}
echo "#SBATCH -e ${path_result}std_error_chicago.txt" >> ${path_result}chicago_${run}
echo "source /beegfs/home/alaverre/.bashrc" >> ${path_result}chicago_${run}
echo "Rscript ~/Tools/CHICAGO/chicagoTools/runChicago.R --export-format interBed,washU_text --design-dir ${path_data} ${files} ${name} -o ${path_result}/" >> ${path_result}chicago_${run}
sbatch -p normal --time=24:00:00 -c 4 --mem=80G ${path_result}chicago_${run}
