#!/bin/bash

sp1=$1
sp2=$2

pathResult="/mnt/result/genome_alignment/${sp1}2other/"
pathScript="/mnt/script/"

for file in `ls ${pathResult}x* | xargs -n 1 basename`
do
python ${pathScript}run_tba.py ${file} ${sp1} ${sp2} &
done

