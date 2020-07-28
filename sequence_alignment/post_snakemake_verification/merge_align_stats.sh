#!/bin/bash

sp_origin=$1
sp_target=$2
data=$3

path=/beegfs/data/alaverre/Regulatory_landscape/test/result/${sp_origin}2other/${sp_origin}2${sp_target}/
path_align=${path}/pecan_alignments

if test -f "${path_align}/AlignmentStatistics_Excluding_all_Exons_${sp_origin}2${sp_target}_${data}.txt"
then
echo "${sp_origin}2${sp_target} ${data} : Merging stats already done !"
sort -u ${path_align}/AlignmentStatistics_Excluding_all_Exons_${sp_origin}2${sp_target}_${data}.txt > ${path_align}/AlignmentStatistics_Excluding_all_Exons_${sp_origin}2${sp_target}_${data}.txt2
grep -v "TotalUngappedLength" ${path_align}/AlignmentStatistics_Excluding_all_Exons_${sp_origin}2${sp_target}_${data}.txt2 > ${path_align}/AlignmentStatistics_Excluding_all_Exons_${sp_origin}2${sp_target}_${data}.txt
sed -i "1i ID.${sp_origin}\tID.${sp_target}\tNbExcludedBases1\tNbExcludedBases2\tTotalUngappedLength\tTotalIdenticalLength\tFilteredUngappedLength\tFilteredIdenticalLength\tTotalAlignmentLength\tFilteredAlignmentLength" ${path_align}/AlignmentStatistics_Excluding_all_Exons_${sp_origin}2${sp_target}_${data}.txt
rm ${path_align}/AlignmentStatistics_Excluding_all_Exons_${sp_origin}2${sp_target}_${data}.txt2

if test -f "${path_align}/AlignmentStatistics_Excluding_all_Exons_${data}.txt_ID_part101" #AlignmentStatistics_Excluding_all_Exons_${data}.txt_ID_part101" or AlignmentStatistics_Excluding_all_Exons_${sp_origin}2${sp_target}_${data}_ID_part100
then
echo "But new stats detected ! So, merging again !"
cat ${path_align}/AlignmentStatistics_Excluding_all_Exons_${data}.txt_ID_part* ${path_align}/AlignmentStatistics_Excluding_all_Exons_${sp_origin}2${sp_target}_${data}.txt > ${path_align}/AlignmentStatistics_Excluding_all_Exons_${sp_origin}2${sp_target}_${data}.txt2
#_${data}.txt_ID_part* or _${sp_origin}2${sp_target}_${data}_ID_part*
grep -v "TotalUngappedLength" ${path_align}/AlignmentStatistics_Excluding_all_Exons_${sp_origin}2${sp_target}_${data}.txt2 | sort -u > ${path_align}/AlignmentStatistics_Excluding_all_Exons_${sp_origin}2${sp_target}_${data}.txt
sed -i "1i ID.${sp_origin}\tID.${sp_target}\tNbExcludedBases1\tNbExcludedBases2\tTotalUngappedLength\tTotalIdenticalLength\tFilteredUngappedLength\tFilteredIdenticalLength\tTotalAlignmentLength\tFilteredAlignmentLength" ${path_align}/AlignmentStatistics_Excluding_all_Exons_${sp_origin}2${sp_target}_${data}.txt
rm ${path_align}/AlignmentStatistics_Excluding_all_Exons_${sp_origin}2${sp_target}_${data}.txt2 ${path_align}/AlignmentStatistics_Excluding_all_Exons_${data}.txt_ID_part* ${path_align}/AlignmentStatistics_Excluding_all_Exons_${sp_origin}2${sp_target}_${data}_ID_part*
fi
else

echo "${sp_origin}2${sp_target} ${data} : Merging..."
cat ${path_align}/AlignmentStatistics_Excluding_all_Exons_${data}.txt_ID_part* > ${path_align}/AlignmentStatistics_Excluding_all_Exons_${sp_origin}2${sp_target}_${data}.txt 
#_${data}.txt_ID_part* or _${sp_origin}2${sp_target}_${data}_ID_part*
grep -v "TotalUngappedLength" ${path_align}/AlignmentStatistics_Excluding_all_Exons_${sp_origin}2${sp_target}_${data}.txt > ${path_align}/AlignmentStatistics_Excluding_all_Exons_${data}.txt2
sed -i "1i ID.${sp_origin}\tID.${sp_target}\tNbExcludedBases1\tNbExcludedBases2\tTotalUngappedLength\tTotalIdenticalLength\tFilteredUngappedLength\tFilteredIdenticalLength\tTotalAlignmentLength\tFilteredAlignmentLength" ${path_align}/AlignmentStatistics_Excluding_all_Exons_${data}.txt2
mv ${path_align}/AlignmentStatistics_Excluding_all_Exons_${data}.txt2 ${path_align}/AlignmentStatistics_Excluding_all_Exons_${sp_origin}2${sp_target}_${data}.txt
rm ${path_align}/AlignmentStatistics_Excluding_all_Exons_${data}.txt_ID_part*
rm ${path_align}/AlignmentStatistics_Excluding_all_Exons_${sp_origin}2${sp_target}_${data}_ID_part*
fi

nb_align=$(wc -l ${path_align}/AlignmentStatistics_Excluding_all_Exons_${sp_origin}2${sp_target}_${data}.txt | awk '{print $1}')
nb_align=$(($nb_align - 1))
nb_total=$(wc -l ${path}/${sp_origin}2${sp_target}_${data}_ID.bed_all_files | awk '{print $1}' )

if [ ${nb_total} -eq ${nb_align} ]
then
echo "Total pairs = Aligned pairs. All done !"
else
echo "Total pairs = ${nb_total} and Aligned pairs = ${nb_align}"
fi

echo ""
