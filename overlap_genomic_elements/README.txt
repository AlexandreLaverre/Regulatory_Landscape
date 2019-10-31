## Overlap interest_file to reference file
./overlap.py ${specie} ${reference_file} ${interest_file} ${output_file}


# list of all overlap done in both mouse and human :
	- merged contacted fragments with : 
		- enhancers_coordinate (CAGE, ENCODE, GRO_seq, RoadMap)
		- exons_coordinate (coding or non_coding genes)
		- TSS or genes coordinates
		- repeats elements
		- phastcons elements

	- lifted merged contacted fragments with same interest files as before but in target specie
