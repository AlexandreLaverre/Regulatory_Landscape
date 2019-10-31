#!/bin/bash

export genome=$1
export i=$2
export cluster=$3
export nthreads=$4

####################################################################################

if [ ${cluster} = "in2p3" ]; then
    export path=/sps/biometr/necsulea/RegulatoryLandscapes
    export pathTools=/sps/biometr/necsulea/Tools
fi

if [ ${cluster} = "cloud" ]; then
    export path=/mnt/RegulatoryLandscapes
    export pathTools=/mnt/Tools
fi

export pathResults=${path}/results/mutual_information_network/${genome}
export pathAracne=${pathTools}/ARACNe-AP/
export pathScripts=${path}/scripts/mutual_information_network

####################################################################################

export regulators=EnhancerList.txt

####################################################################################

if [ -e ${pathResults}/aracne_replicates ]; then
    echo "path out already there"
else
    mkdir ${pathResults}/aracne_replicates 
fi

####################################################################################

if [ ${cluster} = "in2p3" ]; then
    echo "#!/bin/bash" >  ${pathScripts}/bsub_script_aracne
    
    echo "java -Xmx20G -jar ${pathAracne}/dist/aracne.jar -e ${pathResults}/TPM.txt -o ${pathResults}/aracne_replicates --tfs ${pathResults}/${regulators} --pvalue 1E-8 --threads ${nthreads} --seed ${i}" >> ${pathScripts}/bsub_script_aracne
    
    qsub -q mc_huge -l s_rss=6G,sps=1 -pe multicores ${nthreads} -o ${pathScripts}/std_output_aracne_${genome}_${i}.txt -e ${pathScripts}/std_error_aracne_${genome}_${i}.txt ${pathScripts}/bsub_script_aracne

fi

####################################################################################

if [ ${cluster} = "cloud" ]; then
    java -Xmx20G -jar ${pathAracne}/dist/aracne.jar -e ${pathResults}/TPM.txt -o ${pathResults}/aracne_replicates --tfs ${pathResults}/${regulators}  --pvalue 1E-8 --threads ${nthreads} --seed ${i}
fi

####################################################################################
