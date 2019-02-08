#!/bin/bash

####################################################################

export path=/pandata/alaverre
export pathScripts=${path}/scripts
export pathdata=${path}/data/liftover/mouse2other/mouse2human/parts_mouse_enh/simul_uni
export pathresults=${path}/data/liftover/mouse2other/mouse2elephant/results/simul_uni
export pathAlignment=${path}/data/genome_alignments/mouse2other
export pathliftOver=/panhome/alaverre/Tools

####################################################################

for file in `ls ${pathdata}/*bed`
do
export prefix=`basename ${file} .bed`
echo "#!/bin/bash" > ${pathScripts}/bsub_script_liftOver
echo "#PBS -o ${pathresults}/std_output_liftOver.txt" >> ${pathScripts}/bsub_script_liftOver
echo "#PBS -e ${pathresults}/std_error_liftOver.txt" >> ${pathScripts}/bsub_script_liftOver
echo "source /panhome/alaverre/.bashrc" >> ${pathScripts}/bsub_script_liftOver ###PATH mis a jour avec ts programmes
echo "${pathliftOver}/liftOver -minMatch=0.4 ${file} ${pathAlignment}/mm10ToLoxAfr3.over.chain ${pathresults}/${prefix}.lifted.bed ${pathresults}/${prefix}.unmapped " >> ${pathScripts}/bsub_script_liftOver

qsub -q q1hour -l nodes=1:ppn=1,mem=2gb ${pathScripts}/bsub_script_liftOver

done



