###################################################################### Pipeline Sequences Conservation ########################################

##### Lift interest genome coordinates
# 1 - Download map.chain file from http://hgdownload.cse.ucsc.edu/goldenpath/${origin_genome}/liftOver/
# 2 - Run LiftOver 
./liftOver -minMatch=0.1 ${interest_coordinates}_${origin_sp} ../ref/${origin_genome}to${target_genome}.over.chain.gz ${origin_sp}2${target_sp}_${interest_coordinates}.bed ${origin_sp}2${target_sp}_${interest_coordinates}_unMapped

##### Extract sequences of interest coordinates 
# 1 - Download soft masked fasta from Ensembl
# 2 - Extract sequences of interest genome coordinates
./extract_frag_seq.py ${origin_sp} ${target_sp} ${interest_coordinates} # if not lifted coordinates origin_sp == target_sp


##### Align origin coordinates and lifted coordinates with PECAN
# 1 - Make file : list of pairs fragments ID
cut -f 4 ${origin_sp}2${target_sp}_${interest_coordinates}.bed > ID_origin
sed -i 's/$/:+/g' ID_origin
cut -f 1,2,3,6 ${origin_sp}2${target_sp}_${interest_coordinates}.bed > ID_target
sed -i 's/\t/:/g' ID_target
paste ID_origin ID_target > ${origin_sp}2${target_sp}_${interest_coordinates}_list_ID.bed
sed -i 's/ /\t/g' ${origin_sp}2${target_sp}_${interest_coordinates}_list_ID.bed
rm ID_target ID_origin

# 2 - Split it to parralelize
split -l 15000 ${origin_sp}2${target_sp}_${interest_coordinates}_list_ID.bed list/${interest_coordinates}_
for i in `ls list/${interest_coordinates}*`; do sed -i '1i ID.${origin_sp}\tID.{target_sp}' $i; done

# 3 - Run PECAN
run_pecan.sh ${origin_sp} {target_sp} ${interest_coordinates}


#####  Extract exons position
# 1 - Download exons positions from Ensembl : https://www.ensembl.org/biomart/martview with this arguments : 
# Chromosome/scaffold name	Exon region start (bp)	Exon region end (bp)	Exon stable ID	Strand

# 2 - Overlap fragments with exons
./overlap.py ${origin_sp} ${interest_coordinates}_${origin_sp} ${origin_sp}_${exons_file} ${output_file_exons_type}_overlap_${origin_sp}
./overlap.py ${target_sp} ${origin_sp}2${target_sp}_${interest_coordinates}.bed ${target_sp}_${exons_file} ${output_file_exons_type}_overlap_${origin_sp}


##### Extract alignment statistics excluding exons positions
./extract_aln_stats_excluding.regions.sh ${origin_sp} {target_sp} ${interest_coordinates} ${aligner} ${exons_type}



###################################################################### Sequence duplication ###################################################################### 

##### Extract sequences of interest coordinates 
# 1 - Download hard masked fasta from Ensembl
# 2 - Extract sequences of interest genome coordinates
./extract_frag_seq.py ${specie} ${specie} ${interest_coordinates} 

# Split it by chromosome to parralelize : output = ${interest_coordinates}_masked_${chr}.fa

# Run BLAT
./run_blat.sh ${specie} 	# Identity treshold used = 95% 

# Extract BLAT stats with specified treshold of matched length (here we used 80%)
./extract_duplication_stats.py ${specie} ${treshold} ${interest_coordinates}_masked_${chr}.fa 
