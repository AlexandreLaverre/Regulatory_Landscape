#!/bin/bash

export sp=$1
export author=$2
export cluster=$3
export nthreads=$4

##########################################################################

export pathHiCup=/beegfs/data/${USER}/Regulatory_landscape/result/HiCup/${sp}
export pathFinalData=/beegfs/data/necsulea/RegulatoryLandscapesManuscript
export pathFragmentCoords=${pathFinalData}/SupplementaryDataset1/${sp}
export pathIRODS=/lbbeZone/home/alaverre/regulatory_landscape/result/HiCup/${sp}/${author}
export pathLog=/beegfs/data/${USER}/Regulatory_landscape/scripts/estimate_mappability/log/ReadsCoverage_${sp}_${author}

if [ ${sp} = "human" ]; then
    export genome="hg38"
fi

if [ ${sp} = "mouse" ]; then
    export genome="mm10"
fi

##########################################################################
## Formatting file: region coordinates 0-based and whitout ^chr from chromosome names

if [ -e ${pathHiCup}/frag_coords_${genome}_formatted.bed ]; then
    echo "fragment coords already done"
else
    awk '{$2 = $2 - 1; print}' OFS='\t' ${pathFragmentCoords}/frag_coords_${genome}.bed | sed '1d' |  sed -e 's/^chr//'> ${pathHiCup}/frag_coords_${genome}_formatted.bed
fi

##########################################################################

mkdir -p ${pathHiCup}/${author}/${sample}
echo "###################"
echo "${sample}"

for replicat in `ils ${pathIRODS}/${sample} | grep _sorted.dedup.bam`
do
	export prefix=`basename ${replicat} _sorted.dedup.bam`
	echo "${prefix}"

	if [ -e ${pathHiCup}/${author}/${sample}/${prefix}_restriction_fragments_coverage.txt ]; then
	    echo "coverage already done!"
	else

	## Get BAM from IRODS
	if [ -e ${pathHiCup}/${author}/${sample}/${replicat} ]; then
	    echo "BAM already there"
	else
	    echo "Get ${prefix} from IRODS..."
	    iget -N 10 ${pathIRODS}/${sample}/${replicat} ${pathHiCup}/${author}/${sample}
	fi

##########################################################################
	echo "Preparation of BAM file and Calculate reads coverage..."
	echo "#!/bin/bash" > ${pathLog}_${prefix}_sub

	if [ ${cluster} = "pbil" ]; then
		echo "#SBATCH --job-name=ReadsCoverage_${sp}_${author}_${prefix}" >>  ${pathLog}_${prefix}_sub
		echo "#SBATCH --partition=normal" >>  ${pathLog}_${prefix}_sub
		echo "#SBATCH --output=${pathLog}_${prefix}_output" >> ${pathLog}_${prefix}_sub
		echo "#SBATCH --error=${pathLog}_${prefix}_error" >>  ${pathLog}_${prefix}_sub
		echo "#SBATCH --constraint='haswell|skylake|broadwell'" >> ${pathLog}_${prefix}_sub
		echo "#SBATCH --cpus-per-task=${nthreads}" >> ${pathLog}_${prefix}_sub
		echo "#SBATCH --time=10:00:00" >>  ${pathLog}_${prefix}_sub
		echo "#SBATCH --mem=20G" >> ${pathLog}_${prefix}_sub
	fi

##########################################################################
	
	## Sort by position and Index BAM file
	echo "samtools sort -m 2G -@ ${nthreads} -O bam -o ${pathHiCup}/${author}/${sample}/${prefix}_sorted_by_pos.bam ${pathHiCup}/${author}/${sample}/${replicat}" >>  ${pathLog}_${prefix}_sub
	echo "samtools index -@ ${nthreads} ${pathHiCup}/${author}/${sample}/${prefix}_sorted_by_pos.bam" >>  ${pathLog}_${prefix}_sub

	echo "rm ${pathHiCup}/${author}/${sample}/${replicat}" >>  ${pathLog}_${prefix}_sub

	## Calculate coverage
	echo "bedtools coverage -a ${pathHiCup}/frag_coords_${genome}_formatted.bed -b ${pathHiCup}/${author}/${sample}/${prefix}_sorted_by_pos.bam -split -sorted -counts > ${pathHiCup}/${author}/${sample}/${prefix}_restriction_fragments_coverage.txt" >>  ${pathLog}_${prefix}_sub

	#samtools bedcov ${pathHiCup}/frag_coords_${genome}_formatted.bed ${pathHiCup}/${sample}/${prefix}_sorted_by_pos.bam > ${pathHiCup}/${sample}/${prefix}_frag_coverage.txt

	echo "rm ${pathHiCup}/${author}/${sample}/${prefix}*.bam ${pathHiCup}/${author}/${sample}/${prefix}*.bai" >>  ${pathLog}_${prefix}_sub

##########################################################################

	## Running job 
	if [ ${cluster} = "pbil" ]; then
		sbatch ${pathLog}_${prefix}_sub
	else
		bash ${pathLog}_${prefix}_sub
	fi

	## Waiting for output
	until [ -f ${pathHiCup}/${author}/${sample}/${prefix}_restriction_fragments_coverage.txt ]
	do
	     sleep 2m
	done
	echo "${prefix} coverage done!"

	fi
done

##########################################################################

