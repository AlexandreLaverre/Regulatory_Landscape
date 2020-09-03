###################################################################### Sequence duplication ###################################################################### 

##### Extract sequences of interest coordinates 
# 1 - Download hard masked fasta from Ensembl
# 2 - Extract sequences of interest genome coordinates
./extract_rm_seq_dupli.py ${specie}

# Run BLAT
./run_blat.sh ${specie} 	# Identity treshold used = 95% 

# Extract BLAT stats with specified treshold of matched length (here we used 80%)
./extract_duplication_stats.py ${specie} ${treshold} ${interest_coordinates}_masked_${chr}.fa 
