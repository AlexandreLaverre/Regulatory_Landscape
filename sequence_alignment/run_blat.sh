#!/bin/bash
sp=$1
if [ ${sp} = "human" ]; then genome="Homo_sapiens.GRCh38.dna_rm.primary_assembly.fa"; else genome="Mus_musculus.GRCm38.dna_rm.primary_assembly.fa"; fi

path_data="/beegfs/data/alaverre/data/genome_ref/"
path_result="/beegfs/data/alaverre/result/genome_alignment/${sp}2other/duplication_rate/"

for file in `ls ${path_result}${sp}_restriction_fragments_mask/ | xargs -n 1 basename`
do
echo "#!/bin/bash" > ${path_result}running/BLAT_${file}
echo "#SBATCH -o ${path_result}running/std_output_BLAT${file}.txt" >> ${path_result}running/BLAT_${file}
echo "#SBATCH -e ${path_result}running/std_error_BLAT${file}.txt" >> ${path_result}running/BLAT_${file}
echo "source /beegfs/home/alaverre/.bashrc" >> ${path_result}running/BLAT_${file}
echo "blat -out=blast9 -minIdentity=95 ${path_data}${genome} ${path_result}${sp}_restriction_fragments_mask/${file} ${path_result}${file}_output.psl" >> ${path_result}running/BLAT_${file}
sbatch -p normal --time=40:00:00 -c 1 --mem=4G ${path_result}running/BLAT_${file}
done
