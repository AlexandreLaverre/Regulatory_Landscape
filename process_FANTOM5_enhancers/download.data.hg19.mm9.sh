#!/bin/bash

export cloud=$1

##################################################################

if [ ${cloud} = "in2p3" ]; then
    export path=/sps/biometr/necsulea/RegulatoryLandscapes
fi

export pathData=${path}/data/FANTOM5

## downloaded on October 23rd 2019

##################################################################

cd ${pathData}/hg19

wget http://fantom.gsc.riken.jp/5/datafiles/latest/extra/Enhancers/00_readme.txt
wget http://fantom.gsc.riken.jp/5/datafiles/latest/extra/Enhancers/human_permissive_enhancers_phase_1_and_2.bed.gz
wget http://fantom.gsc.riken.jp/5/datafiles/latest/extra/Enhancers/human_permissive_enhancers_phase_1_and_2_expression_tpm_matrix.txt.gz
wget http://fantom.gsc.riken.jp/5/datafiles/latest/extra/Enhancers/Human.sample_name2library_id.txt

wget http://fantom.gsc.riken.jp/5/datafiles/latest/extra/CAGE_peaks/hg19.cage_peak_phase1and2combined_ann.txt.gz
wget http://fantom.gsc.riken.jp/5/datafiles/latest/extra/CAGE_peaks/hg19.cage_peak_phase1and2combined_coord.bed.gz
wget http://fantom.gsc.riken.jp/5/datafiles/latest/extra/CAGE_peaks/hg19.cage_peak_phase1and2combined_counts.osc.txt.gz
wget http://fantom.gsc.riken.jp/5/datafiles/latest/extra/CAGE_peaks/hg19.cage_peak_phase1and2combined_counts_ann.osc.txt.gz
wget http://fantom.gsc.riken.jp/5/datafiles/latest/extra/CAGE_peaks/hg19.cage_peak_phase1and2combined_rel_expr.txt.gz
wget http://fantom.gsc.riken.jp/5/datafiles/latest/extra/CAGE_peaks/hg19.cage_peak_phase1and2combined_tpm.osc.txt.gz
wget http://fantom.gsc.riken.jp/5/datafiles/latest/extra/CAGE_peaks/hg19.cage_peak_phase1and2combined_tpm_ann.osc.txt.gz

##################################################################

cd ${pathData}/mm9

wget http://fantom.gsc.riken.jp/5/datafiles/latest/extra/Enhancers/00_readme.txt

wget http://fantom.gsc.riken.jp/5/datafiles/latest/extra/Enhancers/Mouse.sample_name2library_id.txt
wget http://fantom.gsc.riken.jp/5/datafiles/latest/extra/Enhancers/mouse_permissive_enhancers_phase_1_and_2.bed.gz
wget http://fantom.gsc.riken.jp/5/datafiles/latest/extra/Enhancers/mouse_permissive_enhancers_phase_1_and_2_expression_count_matrix.txt.gz
wget http://fantom.gsc.riken.jp/5/datafiles/latest/extra/Enhancers/mouse_permissive_enhancers_phase_1_and_2_expression_tpm_matrix.txt.gz

wget http://fantom.gsc.riken.jp/5/datafiles/latest/extra/CAGE_peaks/mm9.cage_peak_phase1and2combined_ann.txt.gz
wget http://fantom.gsc.riken.jp/5/datafiles/latest/extra/CAGE_peaks/mm9.cage_peak_phase1and2combined_coord.bed.gz
wget http://fantom.gsc.riken.jp/5/datafiles/latest/extra/CAGE_peaks/mm9.cage_peak_phase1and2combined_counts.osc.txt.gz
wget http://fantom.gsc.riken.jp/5/datafiles/latest/extra/CAGE_peaks/mm9.cage_peak_phase1and2combined_counts_ann.osc.txt.gz
wget http://fantom.gsc.riken.jp/5/datafiles/latest/extra/CAGE_peaks/mm9.cage_peak_phase1and2combined_rel_expr.txt.gz
wget http://fantom.gsc.riken.jp/5/datafiles/latest/extra/CAGE_peaks/mm9.cage_peak_phase1and2combined_tpm.osc.txt.gz
wget http://fantom.gsc.riken.jp/5/datafiles/latest/extra/CAGE_peaks/mm9.cage_peak_phase1and2combined_tpm_ann.osc.txt.gz


##################################################################
