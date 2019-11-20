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

wget -r -np -nd http://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_latest/extra/enhancer

wget -r -np -nd http://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_latest/extra/CAGE_peaks_annotation

wget -r -np -nd http://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_latest/extra/CAGE_peaks_expression

##################################################################

cd ${pathData}/mm10

wget -r -np -nd http://fantom.gsc.riken.jp/5/datafiles/reprocessed/mm10_latest/extra/enhancer

wget -r -np -nd http://fantom.gsc.riken.jp/5/datafiles/reprocessed/mm10_latest/extra/CAGE_peaks_annotation

wget -r -np -nd http://fantom.gsc.riken.jp/5/datafiles/reprocessed/mm10_latest/extra/CAGE_peaks_expression

##################################################################
