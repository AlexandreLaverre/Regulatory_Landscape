#!/bin/bash

export sp1=$1
export sp2=$2

####################################################################

export path=/beegfs/data/alaverre

export pathResults=${path}/result/genome_alignment/mouse2other
export pathScripts=${path}/script 

####################################################################

## pathPairs - tab separated
## ID.Mouse  ID.Rat

## TBA output files have to be called
## ${idmouse}.${idrat}.maf.gz
## they have to be in ${dirTBA}
        
perl ${pathScripts}/extract.aln.stats.pl --species1=${sp1} --species2=${sp2} --pathPairs=${pathResults}/${sp1}2${sp2}_ID.bed   --dirTBA=${pathResults}/tba_alignments/${sp1}2${sp2}/ --minAlignmentLength=10 --pathOutput=${pathResults}/AlignmentStatistics_TBA_${sp1}2${sp2}.txt 

####################################################################



