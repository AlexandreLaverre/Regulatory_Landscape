#!/bin/bash

################################################################################################################################################

ref_sp=$1
target_sp=$2
data=$3

pathScript=/beegfs/data/alaverre/Regulatory_landscape/test/script/
path=/beegfs/data/alaverre/Regulatory_landscape/test/result/${ref_sp}2other/${ref_sp}2${target_sp}
path_result=${path}/pecan_alignments/
path_align=${path_result}/${ref_sp}2${target_sp}_${data}/

################################################################################################################################################

echo "######################################################"
echo "###### Verifying ${ref_sp} to ${target_sp} ${data} alignments ########"
echo "######################################################"

######################### In case of previous check #########################
if test -f "${path}/${ref_sp}2${target_sp}_${data}_ID.bed_all_files"
  then

  if test -f "${path_result}/AlignmentStatistics_Excluding_all_Exons_${ref_sp}2${target_sp}_${data}.txt" # Alignment stats already done
    then
    echo "Alignment stats already done"
    tail -n +2 ${path_result}/AlignmentStatistics_Excluding_all_Exons_${ref_sp}2${target_sp}_${data}.txt | cut -f1  > ${path}/${ref_sp}2${target_sp}_${data}_missing_ID
    grep -Fvf ${path}/${ref_sp}2${target_sp}_${data}_missing_ID  ${path}/${ref_sp}2${target_sp}_${data}_ID.bed_all_files >  ${path}/${ref_sp}2${target_sp}_${data}_missing_pairs
  else
    ls ${path_align} > ${path}/${ref_sp}2${target_sp}_${data}_missing
    cut -d '.' -f1 ${path}/${ref_sp}2${target_sp}_${data}_missing > ${path}/${ref_sp}2${target_sp}_${data}_missing_ID
    grep -Fvf ${path}/${ref_sp}2${target_sp}_${data}_missing_ID  ${path}/${ref_sp}2${target_sp}_${data}_ID.bed_all_files >  ${path}/${ref_sp}2${target_sp}_${data}_missing_pairs
  fi

  nb_missing=$(wc -l ${path}/${ref_sp}2${target_sp}_${data}_missing_pairs | awk '{print $1}' )
  nb_missing_old=$(wc -l ${path}/${ref_sp}2${target_sp}_${data}_ID.bed | awk '{print $1}' )
  echo "A previous check revealed ${nb_missing_old} missing aligment(s)"
  echo "--> New checking detects ${nb_missing} missing aligment(s) !"

  if [ ${nb_missing} -gt 0 ]
    then
    mv ${path}/${ref_sp}2${target_sp}_${data}_missing_pairs ${path}/${ref_sp}2${target_sp}_${data}_ID.bed
    rm ${path}/${ref_sp}2${target_sp}_${data}_missing ${path}/${ref_sp}2${target_sp}_${data}_missing_ID
    rm ${path}/list/*${data}*
    rm -r ${path}/pecan_alignments/running_${data}/
    mkdir ${path}/pecan_alignments/running_${data}/

    echo "Running alignment pipeline !"
    part=$((${nb_missing} / 5000 ))
    part_int=$((${part%.*} + 1 ))
    screen -Sdm Snakemake_${ref_sp}2${target_sp}_${data}_missing ${pathScript}/run_Snakemake.sh ${ref_sp} ${target_sp} ${data} ${part_int}
  else
    echo "It's all good ! No missing files !"
    rm ${path}/${ref_sp}2${target_sp}_${data}_missing_pairs ${path}/${ref_sp}2${target_sp}_${data}_missing_ID
  fi

######################### Without previous check #########################
else
  if test -f "${path_result}/AlignmentStatistics_Excluding_all_Exons_${ref_sp}2${target_sp}_${data}.txt" # Alignment extract stats already done
    then
    echo "Alignment stats already done"
    tail -n +2 ${path_result}/AlignmentStatistics_Excluding_all_Exons_${ref_sp}2${target_sp}_${data}.txt | cut -f1  > ${path}/${ref_sp}2${target_sp}_${data}_missing_ID
    grep -Fvf ${path}/${ref_sp}2${target_sp}_${data}_missing_ID  ${path}/${ref_sp}2${target_sp}_${data}_ID.bed >  ${path}/${ref_sp}2${target_sp}_${data}_missing_pairs

  else
      if test -d "${path_align}" # Alignment already done
        then
        ls ${path_align} > ${path}/${ref_sp}2${target_sp}_${data}_missing
        cut -d '.' -f1 ${path}/${ref_sp}2${target_sp}_${data}_missing > ${path}/${ref_sp}2${target_sp}_${data}_missing_ID
        grep -Fvf ${path}/${ref_sp}2${target_sp}_${data}_missing_ID  ${path}/${ref_sp}2${target_sp}_${data}_ID.bed >  ${path}/${ref_sp}2${target_sp}_${data}_missing_pairs
      else
        echo "No alignment found -->  Running alignment pipeline !"
        screen -Sdm Snakemake_${ref_sp}2${target_sp}_${data}_missing ${pathScript}/run_Snakemake.sh ${ref_sp} ${target_sp} ${data} 50
        echo ""
        exit
      fi
  fi

  nb_missing=$(wc -l ${path}/${ref_sp}2${target_sp}_${data}_missing_pairs | awk '{print $1}' )

  mv ${path}/${ref_sp}2${target_sp}_${data}_ID.bed ${path}/${ref_sp}2${target_sp}_${data}_ID.bed_all_files
  mv ${path}/${ref_sp}2${target_sp}_${data}_missing_pairs ${path}/${ref_sp}2${target_sp}_${data}_ID.bed
  rm ${path}/${ref_sp}2${target_sp}_${data}_missing ${path}/${ref_sp}2${target_sp}_${data}_missing_ID

  if [ ${nb_missing} -gt 0 ]
    then
    echo "${ref_sp}2${target_sp} in ${data} : Missing ${nb_missing} aligment(s) !"
    rm -r ${path}/pecan_alignments/running_${data}/
    mkdir ${path}/pecan_alignments/running_${data}/

    echo "Running alignment pipeline !"
    part=$((${nb_missing} / 5000 ))
    part_int=$((${part%.*} + 1 ))
    screen -Sdm Snakemake_${ref_sp}2${target_sp}_${data}_missing ${pathScript}/run_Snakemake.sh ${ref_sp} ${target_sp} ${data} ${part_int}
  else
    echo "It's all good ! No missing files !"
    rm ${path}/${ref_sp}2${target_sp}_${data}_missing_pairs ${path}/${ref_sp}2${target_sp}_${data}_missing_ID
  fi
fi

echo ""

### export sp=(mouse cow opossum elephant rat rabbit macaque dog chicken)
### for i in "${sp[@]}"; do ./check_alignments.sh human $i CAGE; done 
