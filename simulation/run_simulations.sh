#!/bin/bash

sp=$1 # ex: human or mouse
data=$2 # ex : bait_all or bait_other

pathResult="/beegfs/data/alaverre/Regulatory_landscape/result/simulations/${sp}_samples/"
pathScript="/beegfs/data/alaverre/Regulatory_landscape/script/"


if [ $sp == "human" ]
then
cells=("TCD4Non" "TCD4MF" "TCD8" "TCD4Act" "TB" "pre_adipo" "PEK_undiff" "PEK_late" "PEK_early" "Neu" "NCD8" "NCD4" "NB" "Mon" "MK" "Mac2" "Mac1" "Mac0" "hNEC" "hESC" "FoeT" "Ery" "EP" "CD34" "cardio" "Bcell")
genome="hg38"
else
cells=("EpiSC" "ESC_18" "ESC" "ESC_NKO" "ESC_wild" "ESd_starved" "ESd_TPO" "FLC" "preadip_4H" "preadip_D0" "preadip_D2" "preB_aged" "preB_young" "TSC")
genome="mm10"
fi

for cell in "${cells[@]}"
do
echo "#!/bin/bash" > ${pathResult}simulations_${data}_${sp}_${cell}
echo "#SBATCH -o ${pathResult}std_output_${data}_${cell}.txt" >> ${pathResult}simulations_${data}_${sp}_${cell}
echo "#SBATCH -e ${pathResult}std_error_${data}_${cell}.txt" >> ${pathResult}simulations_${data}_${sp}_${cell}
echo "source /beegfs/home/alaverre/.bashrc" >> ${pathResult}simulations_${data}_${sp}_${cell}
echo "python ${pathScript}simulations_fragoverbin_${data}_pbil.py ${sp} ${genome} ${cell}" >> ${pathResult}simulations_${data}_${sp}_${cell}
sbatch -p normal --time=30:00:00 --mem=1GB -c 1 ${pathResult}simulations_${data}_${sp}_${cell}
done


