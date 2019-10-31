#!/bin/bash

export genome=$1
export cluster=$2

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

if [ ${cluster} = "in2p3" ]; then

    echo "#!/bin/bash" >  ${pathScripts}/bsub_script_consolidate
    
    echo "java -Xmx20G -jar ${pathAracne}/dist/aracne.jar -o ${pathResults}/aracne_replicates --consolidate --nobonferroni" >>  ${pathScripts}/bsub_script_consolidate

    echo "mv ${pathResults}/aracne_replicates/network.txt  ${pathResults}/aracne_replicates/network_withoutBonferroni.txt"  >>  ${pathScripts}/bsub_script_consolidate

     echo "java -Xmx20G -jar ${pathAracne}/dist/aracne.jar -o ${pathResults}/aracne_replicates --consolidate " >>  ${pathScripts}/bsub_script_consolidate

    echo "mv ${pathResults}/aracne_replicates/network.txt  ${pathResults}/aracne_replicates/network_withBonferroni.txt"  >>  ${pathScripts}/bsub_script_consolidate
    
    qsub -q mc_highmem_huge -l s_rss=30G,sps=1 -o ${pathScripts}/std_output_aracne_consolidate_${genome}.txt -e ${pathScripts}/std_error_aracne_${consolidate}_${genome}.txt ${pathScripts}/bsub_script_consolidate
  
fi

####################################################################################

if [ ${cluster} = "cloud" ]; then
    java -Xmx20G -jar ${pathAracne}/dist/aracne.jar -o ${pathResults}/aracne_replicates --consolidate --nobonferroni
    
    mv ${pathResults}/aracne_replicates/network.txt  ${pathResults}/aracne_replicates/network_withoutBonferroni.txt
    
    java -Xmx20G -jar ${pathAracne}/dist/aracne.jar -o ${pathResults}/aracne_replicates --consolidate
    
    mv ${pathResults}/aracne_replicates/network.txt  ${pathResults}/aracne_replicates/network_withBonferroni.txt
fi

####################################################################################
