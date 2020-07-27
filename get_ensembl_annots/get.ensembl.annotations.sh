#######################################################################

export species=$1
export release=$2

#######################################################################

export pathResults=../../../RegulatoryLandscapesManuscript/SupplementaryDataset3

#######################################################################

if [ ${species} = "mouse" ]; then
    export db="mus_musculus_core_${release}_38"
fi

if [ ${species} = "rat" ]; then
    export db="rattus_norvegicus_core_${release}_6"
fi

if [ ${species} = "human" ]; then
    export db="homo_sapiens_core_${release}_38"
fi

if [ ${species} = "macaque" ]&&[ ${release} = "99" ]; then
    export db="macaca_mulatta_core_${release}_10"
fi

if [ ${species} = "chicken" ]; then
    export db="gallus_gallus_core_${release}_5"
fi

if [ ${species} = "cow" ]; then
    export db="bos_taurus_core_${release}_31"
fi

if [ ${species} = "dog" ]; then
    export db="canis_familiaris_core_${release}_31"
fi

if [ ${species} = "elephant" ]; then
    export db="loxodonta_africana_core_${release}_3"
fi

if [ ${species} = "opossum" ]; then
    export db="monodelphis_domestica_core_${release}_5"
fi

if [ ${species} = "rabbit" ]; then
    export db="oryctolagus_cuniculus_core_${release}_2"
fi

#######################################################################

## transcript info

echo "use $db; " > mysql.script.sh
echo "select gene.stable_id, transcript.stable_id, transcript.biotype, seq_region.name, transcript.seq_region_start, transcript.seq_region_end, transcript.seq_region_strand from transcript, gene, seq_region where transcript.gene_id=gene.gene_id and transcript.seq_region_id=seq_region.seq_region_id; " >> mysql.script.sh

echo GeneID$'\t'TranscriptID$'\t'TranscriptBiotype$'\t'Chr$'\t'Start$'\t'End$'\t'Strand >${pathResults}/transcripts/${species}_transcripts_Ensembl${release}.txt
mysql -h ensembldb.ensembl.org -P 5306 -u anonymous < mysql.script.sh | sed '1d' >> ${pathResults}/transcripts/${species}_transcripts_Ensembl${release}.txt

#######################################################################

## gene info 

echo "use $db; " > mysql.script.sh
echo "select gene.stable_id, gene.biotype, seq_region.name, gene.seq_region_start, gene.seq_region_end, gene.seq_region_strand from gene, seq_region where gene.seq_region_id=seq_region.seq_region_id; " >> mysql.script.sh

echo GeneID$'\t'GeneBiotype$'\t'Chr$'\t'Start$'\t'End$'\t'Strand>${pathResults}/genes/${species}_genes_Ensembl${release}.txt
mysql -h ensembldb.ensembl.org -P 5306 -u anonymous < mysql.script.sh | sed '1d' >> ${pathResults}/genes/${species}_genes_Ensembl${release}.txt

#######################################################################

rm mysql.script.sh

#######################################################################
