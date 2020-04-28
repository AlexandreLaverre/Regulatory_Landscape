#!/bin/bash
sp=$1
seq=$2

if [ ${sp} = "human" ]; then genome="hg38_masked.fa"; else genome="mm10_masked.fa"; fi

path_data="/beegfs/data/alaverre/Regulatory_landscape/test/data/genomes_rm/"
path_result="/beegfs/data/alaverre/Regulatory_landscape/test/result/${sp}2other/duplication_rate/"

for file in `ls ${path_result}${sp}_restriction_fragments_mask/*${seq}* | xargs -n 1 basename`
do
echo "#!/bin/bash" > ${path_result}running/BLAT_${file}
echo "#SBATCH -o ${path_result}running/std_output_BLAT${file}.txt" >> ${path_result}running/BLAT_${file}
echo "#SBATCH -e ${path_result}running/std_error_BLAT${file}.txt" >> ${path_result}running/BLAT_${file}
echo "source /beegfs/home/alaverre/.bashrc" >> ${path_result}running/BLAT_${file}
echo "blat -out=psl -minIdentity=95 ${path_data}${genome} ${path_result}${sp}_restriction_fragments_mask/${file} ${path_result}/output.psl/${file}_output.psl" >> ${path_result}running/BLAT_${file}
sbatch -p normal --time=48:00:00 -c 1 --mem=4G ${path_result}running/BLAT_${file}
done
