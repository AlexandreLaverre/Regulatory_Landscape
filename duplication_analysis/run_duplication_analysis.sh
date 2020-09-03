#!/bin/bash

sp=$1 # ie: human mouse
seq=$2 # ie: restriction_fragments ENCODE ...
treshold=$3 # ie : 0.8

path="/beegfs/data/alaverre/Regulatory_landscape/test/"
pathScript=path+"script/"
pathResult=path+"result/${sp}2other/duplication_rate/"

######  Extracting masked sequences ##### 

${pathScript}/extract_rm_seq_dupli.py ${sp} ${seq}

######  Running BLAT ##### 

${pathScript}/run_blat.sh ${sp} ${seq}

###### Extract statistics ###### 

for i in `ls ${pathResult}/${sp}_restriction_fragments_mask/`
do
${pathScript}/extract_duplication_stats.py ${sp} ${treshold} ${i}
&
done

###### Merging results ###### 

cat ${pathResult}/*${seq}*txt_GC > ${pathResult}/${sp}_${seq}_BLAT_summary_${treshold}.txt
sed -i "1i frag_ID\tlength\tno_N\tnb_GC\tnb_match\tmatch" ${pathResult}/${sp}_${seq}_BLAT_summary_${treshold}.txt
#rm ${pathResult}/*${seq}*txt_GC

