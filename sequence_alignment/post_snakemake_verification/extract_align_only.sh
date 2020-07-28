#!/bin/bash

sp_origin=$1
sp_target=$2
data=$3
part=$4

echo "This is part : ${part}"

path_pairs=/beegfs/data/alaverre/Regulatory_landscape/test/result/${sp_origin}2other/${sp_origin}2${sp_target}/list
path_exon=/beegfs/data/alaverre/Regulatory_landscape/test/data/${sp_origin}2other/overlap_frag_exon
path_align=/beegfs/data/alaverre/Regulatory_landscape/test/result/${sp_origin}2other/${sp_origin}2${sp_target}/pecan_alignments
path_output=/beegfs/data/alaverre/Regulatory_landscape/test/result/${sp_origin}2other/${sp_origin}2${sp_target}/pecan_alignments


perl ./extract.aln.stats.excluding.regions.pl --species1=${sp_origin} --species2=${sp_target} --pathPairs=${path_pairs}/${part} --pathExcludedRegionsSpecies1=${path_exon}/${sp_origin}_${data}_overlap_all_exons.txt --pathExcludedRegionsSpecies2=${path_exon}/${sp_origin}2${sp_target}_${data}_overlap_all_exons.txt --dirAlignments=${path_align}/${sp_origin}2${sp_target}_${data} --suffixAlignments=.mfa.gz --format=MFA --minAlignmentLength=10 --pathOutput=${path_output}/AlignmentStatistics_Excluding_all_Exons_${sp_origin}2${sp_target}_${part}
