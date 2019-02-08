#!/bin/bash
specie=$1
genome=$2
author_ctype=$3

path_data="/pandata/alaverre/data/CHICAGO_files/${specie}/${genome}/"
path_result="/pandata/alaverre/result/CHICAGO/${specie}/${author_ctype}/"

sample=$(ls ${path_result}*.chinput | sed 's/\t/;/g')
files=$(echo ${sample} | sed 's/ /,/g')

Rscript ~/Tools/CHICAGO/chicagoTools/runChicago.R --export-format interBed,washU_text --design-dir ${path_data} ${files} ESC -o ${path_result}
