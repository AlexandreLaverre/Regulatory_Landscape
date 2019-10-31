#!/bin/sh

sp1=$1  		# ex : human
sp2=$2			# ex : mouse
file=$3			# ex : CAGE, merged
aligner=$4 		# TBA or PECAN
overlap=$5		# nocoding, coding or all

if [ ${aligner} = "TBA" ]
then 
format="MAF"
suffix=".maf.gz"
else
format="MFA"
suffix=".mfa.gz"
fi

####################################################################

path=/beegfs/data/alaverre/Regulatory_landscape
pathResult=${path}/result/genome_alignment/${sp1}2other/
pathScript=${path}/script/

####################################################################

## pathPairs - tab separated
## ID.Mouse  ID.Rat

## TBA output files have to be called
## ${idmouse}.${idrat}.maf.gz
## they have to be in ${dirTBA}

## pathExcludedRegionsSpecies1 and pathExcludedRegionsSpecies2
## regions to be excluded: current format is ok, e.g. 1:65419:65433,1:65520:65573,1:69037:70108
## has to have a header: ID chr start end overlap_ID (tab-separated)

for list in `ls ${pathResult}${sp1}2${sp2}/list/${file}* | xargs -n 1 basename` # list of pairs ${file}*
do
echo "#!/bin/sh" > ${pathResult}${sp1}2${sp2}/pecan_alignments/extract.aln.stats_${sp1}2${sp2}_${list} 
echo "#SBATCH -o ${pathResult}${sp1}2${sp2}/pecan_alignments/std_output_${list}.txt" >> ${pathResult}${sp1}2${sp2}/pecan_alignments/extract.aln.stats_${sp1}2${sp2}_${list} 
echo "#SBATCH -e ${pathResult}${sp1}2${sp2}/pecan_alignments/std_error_${list}.txt" >> ${pathResult}${sp1}2${sp2}/pecan_alignments/extract.aln.stats_${sp1}2${sp2}_${list} 	
echo "source /beegfs/home/alaverre/.bashrc" >> ${pathResult}${sp1}2${sp2}/pecan_alignments/extract.aln.stats_${sp1}2${sp2}_${list} 

echo "perl ${pathScripts}extract.aln.stats.excluding.regions.pl --species1=${sp1} --species2=${sp2} --pathPairs=${pathResult}${sp1}2${sp2}/list/${list} --pathExcludedRegionsSpecies1=${pathResult}overlap_frag_exon/${sp1}_${file}_overlap_${overlap}_exons.txt --pathExcludedRegionsSpecies2=${pathResult}overlap_frag_exon/${sp1}2${sp2}_${file}_overlap_${overlap}_exons.txt  --dirAlignments=${pathResult}${sp1}2${sp2}/pecan_alignments/${sp1}2${sp2}_${file} --suffixAlignments=${suffix} --format=${format} --minAlignmentLength=10  --pathOutput=${pathResult}${sp1}2${sp2}/pecan_alignments/AlignmentStatistics_${aligner}_Excluding_${overlap}_Exons.txt_${list}" >> ${pathResult}${sp1}2${sp2}/pecan_alignments/extract.aln.stats_${sp1}2${sp2}_${list} 

sbatch -p normal --time=1:00:00 --mem=1GB -c 1 ${pathResult}${sp1}2${sp2}/pecan_alignments/extract.aln.stats_${sp1}2${sp2}_${list}

done
####################################################################
