############################################### Pipeline : reads to significant interactions in PC-HiC ###############################################

#### Download SRA files ####
for file in `cat SRA_list.txt`
do
prefetch ${file} 
done

### SRA to FastQ ###
for file in `ls ${path_data}/reads/sra/*.sra | xargs -n 1 basename`
do
echo "#!/bin/bash" > bsub_fastq_${file}
echo "#PBS -o ${path_data}/reads_PC_HIC/${specie}/${data}/std_output.txt" >> bsub_fastq_${file}
echo "#PBS -e ${path_data}/reads_PC_HIC/${specie}/${data}/std_error.txt" >> bsub_fastq_${file}
echo "source .bashrc" >> bsub_fastq_${file}
echo "parallel-fastq-dump --split-files --gzip --threads 10 --sra-id  ${path_data}/reads/sra/${file} --outdir  ${path_data}/reads_PC_HIC/${specie}/${data}/${file} --tmpdir ${path_data}/reads_PC_HIC/${specie}/${data}/${file}" >> bsub_fastq_${file}
mkdir -p ${path_data}/reads_PC_HIC/${specie}/${data}/${file}
qsub -q q1day -l nodes=1:ppn=10,mem=1gb bsub_fastq_${file}
done

### Checking : Count reads number in fastq ####
for i in `ls *fastq.gz`
do 
echo ${i}
zcat ${i} | echo $((`wc -l`/4))
done


#### Create indexed genome #### (quick)
bowtie2-build ${species}.dna.primary_assembly.fa ${genome_name}


#### HiCUP Pipe #### (quick)
# 1 - Genome digest from restriction enzyme
hicup_digester --genome ${genome_name} --re1 A^AGCTT,HindIII *fa 

# 2 - HiCUP 
run_HiCup.sh ${species} ${genome} ${data}

# 3 - Merge technical replicat BAM for each biological replicat
run_merge_dedup ${specie} ${dataset} ${cell_rep}


#### CHICAGO Pipe ######
## CHICAGO tools : preparing files ##
# 1 - baitmap + rmap #
create_baitmap_rmap.pl digest_genome.txt bait_position.txt # bait_position extract from supplementary data in Schoenfelder

# 2 - design files #
python ~/Tools/CHICAGO/chicagoTools/makeDesignFiles.py --rmapfile=${path_data}/CHICAGO_files/${specie}/${genome}/auto/Digest_${genome}_HindIII_None.txt.rmap --baitmapfile=${path_data}/CHICAGO_files/${specie}/${genome}/auto/Digest_${genome}_HindIII.txt.baitmap --outfilePrefix=${path_data}/CHICAGO_files/{specie}/${genome}/auto/${genome}

# 3 - processing BAM files into Chicago input files #
bam2chicago.sh HiCUP/${cell_rep}.hicup.bam Digest_${genome}_HindIII.txt.baitmap Digest_${genome}_HindIII_None.txt.rmap ${cell_rep}

# 4 - CHICAGO
run_chicago.sh ${sp} ${genome} ${dataset} ${cell_rep}


#### Unique interactions file per specie
run uniq_interactions_file_${sp}.py

### Merging contacted sequences by bait
run Merging_contacts.py



