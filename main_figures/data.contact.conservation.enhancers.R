#######################################################################################

library(data.table)

source("parameters.R")

pathEvolution=paste(pathFinalData, "SupplementaryDataset7", sep="")

load(paste(pathFigures, "RData/data.gene.enhancer.contacts.RData", sep=""))
load(paste(pathFigures, "RData/data.ortho.genes.RData", sep=""))

#######################################################################################

contact.conservation=list()

for(ref in c("human", "mouse")){

  tg=setdiff(c("human", "mouse"), ref)

  contact.conservation[[paste(ref, "2", tg, sep="")]]=list()

  for(enh in enhancer.datasets[[ref]]){
    message(enh)
    
    ## filtered gene enhancer contacts

    filtered.contacts.obs=gene.enhancer.contacts[[ref_sp]][[enh]][["real"]]
    filtered.contacts.sim=gene.enhancer.contacts[[ref_sp]][[enh]][["simulated"]]

    filtered.contacts.obs$id=paste(filtered.contacts.obs$gene, filtered.contacts.obs$enhancer, sep="-")
    filtered.contacts.sim$id=paste(filtered.contacts.sim$gene, filtered.contacts.sim$enhancer, sep="-")
        
    obs=fread(paste(pathEvolution, "/", ref, "/contact_conservation/", enh, "/", ref, "_original2", tg,"_original.txt", sep=""), h=T, stringsAsFactors=F, sep="\t")
    sim=fread(paste(pathEvolution, "/", ref, "/contact_conservation/", enh, "/", ref, "_simulated2", tg,"_simulated.txt", sep=""), h=T, stringsAsFactors=F, sep="\t")

    class(obs)="data.frame"
    class(sim)="data.frame"
    
    print(paste(nrow(obs)," observed contacts before filtering"))
    print(paste(nrow(sim)," simulated contacts before filtering"))

    ## select only previously filtered gene - enhancer pairs

    obs$id=paste(obs$origin_gene, obs$origin_enh, sep="-")
    sim$id=paste(sim$origin_gene, sim$origin_enh, sep="-")

    obs=obs[which(obs$id%in%filtered.contacts.obs$id),]
    sim=sim[which(sim$id%in%filtered.contacts.sim$id),]

    print(paste(nrow(obs)," observed contacts after filtering"))
    print(paste(nrow(sim)," simulated contacts after filtering"))
                
  
    ## filtered enhancers - single BLAT match, within accepted distance range
    obs=obs[which(obs$BLAT_match == 1 & obs$origin_dist >= minDistance & obs$origin_dist <= maxDistance),]
    sim=sim[which(sim$BLAT_match == 1 & sim$origin_dist >= minDistance & sim$origin_dist <= maxDistance),]

    ## take only well conserved enhancers
    align.threshold <- 0.4 

    obs=obs[which(obs$align_score>=align.threshold),]
    sim=sim[which(sim$align_score>=align.threshold),]

    ## take only orthologous genes presents in both species datasets
    
    obs=obs[which(obs$target_data == "TRUE"),]
    sim=sim[which(sim$target_data == "TRUE"),]
        
    ## select previously filtered ortho genes
    
    obs=obs[which(obs$origin_gene%in%ortho[,ref] & obs$target_gene%in%ortho[,tg]),]
    sim=sim[which(sim$origin_gene%in%ortho[,ref] & sim$target_gene%in%ortho[,tg]),]
     
    ## save data
    
    contact.conservation[[paste(ref, "2", tg, sep="")]][[enh]]=list("obs"=obs, "sim"=sim) 
    
  }
}

#######################################################################################

save(list=c("contact.conservation"), file=paste(pathFigures, "RData/data.contact.conservation.enhancers.RData",sep=""))

#######################################################################################

