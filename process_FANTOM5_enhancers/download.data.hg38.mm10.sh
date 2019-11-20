#!/bin/bash

export cloud=$1

##################################################################

if [ ${cloud} = "in2p3" ]; then
    export path=/sps/biometr/necsulea/RegulatoryLandscapes
fi

export pathData=${path}/data/FANTOM5

## downloaded on October 2nd 2019
## corresponds to v7, from 15th March 2019

##################################################################

cd ${pathData}/hg38

## enhancers
wget http://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_latest/extra/enhancer/F5.hg38.enhancers.bed.gz
wget http://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_latest/extra/enhancer/F5.hg38.enhancers.expression.matrix.gz
wget http://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_latest/extra/enhancer/F5.hg38.enhancers.expression.tpm.matrix.gz
wget http://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_latest/extra/enhancer/F5.hg38.enhancers.expression.usage.matrix.gz 

## CAGE annotation
wget http://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_latest/extra/CAGE_peaks_annotation/hg38_liftover+new_CAGE_peaks_phase1and2_annot.txt.gz
wget http://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_latest/extra/CAGE_peaks_annotation/hg38_liftover+new_CAGE_peaks_phase1and2_trans.txt.gz
wget http://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_latest/extra/CAGE_peaks_annotation/human_phase1and2_CAGE_Peak_name.txt.gz

## CAGE peaks expression
wget http://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_latest/extra/CAGE_peaks_expression/hg38_fair+new_CAGE_peaks_phase1and2_ann.txt.gz
wget http://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_latest/extra/CAGE_peaks_expression/hg38_fair+new_CAGE_peaks_phase1and2_counts.osc.txt.gz
wget http://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_latest/extra/CAGE_peaks_expression/hg38_fair+new_CAGE_peaks_phase1and2_counts_ann.osc.txt.gz
wget http://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_latest/extra/CAGE_peaks_expression/hg38_fair+new_CAGE_peaks_phase1and2_tpm.osc.txt.gz
wget http://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_latest/extra/CAGE_peaks_expression/hg38_fair+new_CAGE_peaks_phase1and2_tpm_ann.osc.txt.gz
wget http://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_latest/extra/CAGE_peaks_expression/hg38_liftover+new_CAGE_peaks_phase1and2_ann.txt.gz

##################################################################

cd ${pathData}/mm10

## enhancers
wget http://fantom.gsc.riken.jp/5/datafiles/reprocessed/mm10_latest/extra/enhancer/F5.mm10.enhancers.bed.gz
wget http://fantom.gsc.riken.jp/5/datafiles/reprocessed/mm10_latest/extra/enhancer/F5.mm10.enhancers.expression.matrix.gz
wget http://fantom.gsc.riken.jp/5/datafiles/reprocessed/mm10_latest/extra/enhancer/F5.mm10.enhancers.expression.usage.matrix.gz

## CAGE annotation

wget http://fantom.gsc.riken.jp/5/datafiles/reprocessed/mm10_latest/extra/CAGE_peaks_annotation/mm10_liftover+new_CAGE_peaks_phase1and2_annot.txt.gz
wget http://fantom.gsc.riken.jp/5/datafiles/reprocessed/mm10_latest/extra/CAGE_peaks_annotation/mm10_liftover+new_CAGE_peaks_phase1and2_trans.txt.gz
wget http://fantom.gsc.riken.jp/5/datafiles/reprocessed/mm10_latest/extra/CAGE_peaks_annotation/mouse_phase1and2_CAGE_Peak_name.txt.gz

## CAGE peaks expression
wget http://fantom.gsc.riken.jp/5/datafiles/reprocessed/mm10_latest/extra/CAGE_peaks_expression/mm10_fair+new_CAGE_peaks_phase1and2_ann.txt.gz
wget http://fantom.gsc.riken.jp/5/datafiles/reprocessed/mm10_latest/extra/CAGE_peaks_expression/mm10_fair+new_CAGE_peaks_phase1and2_counts.osc.txt.gz
wget http://fantom.gsc.riken.jp/5/datafiles/reprocessed/mm10_latest/extra/CAGE_peaks_expression/mm10_fair+new_CAGE_peaks_phase1and2_counts_ann.osc.txt.gz
wget http://fantom.gsc.riken.jp/5/datafiles/reprocessed/mm10_latest/extra/CAGE_peaks_expression/mm10_fair+new_CAGE_peaks_phase1and2_tpm.osc.txt.gz
wget http://fantom.gsc.riken.jp/5/datafiles/reprocessed/mm10_latest/extra/CAGE_peaks_expression/mm10_fair+new_CAGE_peaks_phase1and2_tpm_ann.osc.txt.gz
wget http://fantom.gsc.riken.jp/5/datafiles/reprocessed/mm10_latest/extra/CAGE_peaks_expression/mm10_liftover+new_CAGE_peaks_phase1and2_ann.txt.gz

##################################################################
