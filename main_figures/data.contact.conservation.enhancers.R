#######################################################################################

library(data.table)

source("parameters.R")

pathEvolution=paste(pathFinalData, "SupplementaryDataset7", sep="")

load(paste(pathFigures, "RData/data.gene.enhancer.contacts.RData", sep=""))
load(paste(pathFigures, "RData/data.ortho.genes.RData", sep=""))
load(paste(pathFigures, "RData/data.bait.annotation.RData", sep=""))
load(paste(pathFigures, "RData/data.gene.annotation.RData", sep=""))

## alignment score minAlignScore defined in parameters.R

#######################################################################################

unfiltered.contact.conservation=list()
contact.conservation=list()

for(ref in c("human", "mouse")){

  tg=setdiff(c("human", "mouse"), ref)

  annot.tg=gene.annot[[tg]]

  contact.conservation[[paste(ref, "2", tg, sep="")]]=list()

  ## baited genes in target species

  baited.genes.tg=baited.genes[[tg]]

  for(enh in enhancer.datasets[[ref]]){
    message(enh)

    ## load sequence conservation 
    load(paste(pathFigures, "RData/data.sequence.conservation.enhancers.",enh,".",ref,"2", tg,".RData", sep=""))
    
    ## filtered gene enhancer contacts

    filtered.contacts.obs=gene.enhancer.contacts[[ref]][[enh]][["real"]]
    filtered.contacts.sim=gene.enhancer.contacts[[ref]][[enh]][["simulated"]]

    filtered.contacts.obs$id=paste(filtered.contacts.obs$gene, filtered.contacts.obs$enhancer, sep="-")
    filtered.contacts.sim$id=paste(filtered.contacts.sim$gene, filtered.contacts.sim$enhancer, sep="-")
        
    obs=fread(paste(pathEvolution, "/", ref, "/contact_conservation/", enh, "/", ref, "_original2", tg,"_original.txt", sep=""), h=T, stringsAsFactors=F, sep="\t")
    sim=fread(paste(pathEvolution, "/", ref, "/contact_conservation/", enh, "/", ref, "_simulated2", tg,"_simulated.txt", sep=""), h=T, stringsAsFactors=F, sep="\t")

    class(obs)="data.frame"
    class(sim)="data.frame"

    nbgenes.obs=length(unique(obs$origin_gene))
    nbgenes.sim=length(unique(sim$origin_gene))
    
    print(paste(nrow(obs)," observed contacts for",nbgenes.obs,"genes before filtering"))
    print(paste(nrow(sim)," simulated contacts for",nbgenes.sim,"genes before filtering"))

    ## select only previously filtered gene - enhancer pairs

    obs$id=paste(obs$origin_gene, obs$origin_enh, sep="-")
    sim$id=paste(sim$origin_gene, sim$origin_enh, sep="-")

    obs=obs[which(obs$id%in%filtered.contacts.obs$id),]
    sim=sim[which(sim$id%in%filtered.contacts.sim$id),]

    nbgenes.obs=length(unique(obs$origin_gene))
    nbgenes.sim=length(unique(sim$origin_gene))

    print(paste(nrow(obs)," observed contacts for",nbgenes.obs,"genes after basic filtering (keeping previously filtered contacts)"))
    print(paste(nrow(sim)," simulated contacts for",nbgenes.sim,"genes after basic filtering  (keeping previously filtered contacts)"))
    
    ## select previously filtered ortho genes
    
    obs=obs[which(obs$origin_gene%in%ortho[,ref] & obs$target_gene%in%ortho[,tg]),]
    sim=sim[which(sim$origin_gene%in%ortho[,ref] & sim$target_gene%in%ortho[,tg]),]

    nbgenes.obs=length(unique(obs$origin_gene))
    nbgenes.sim=length(unique(sim$origin_gene))

    print(paste(nrow(obs)," observed contacts for",nbgenes.obs,"genes after filtering ortho"))
    print(paste(nrow(sim)," simulated contacts for",nbgenes.sim,"genes after filtering ortho"))

    ## origin distance
    ## we use the minimum distance between baited TSS and enhancer
    obs$origin_dist=obs$origin_min_dist_baitedTSS
    sim$origin_dist=sim$origin_min_dist_baitedTSS
    
    ## save data after the first filtering steps

    unfiltered.contact.conservation[[paste(ref, "2", tg, sep="")]][[enh]]=list("obs"=obs, "sim"=sim) 
    
    ## take only genes baited in both species datasets
    
    obs=obs[which(obs$target_gene %in% baited.genes.tg),]
    sim=sim[which(sim$target_gene %in% baited.genes.tg),]

    nbgenes.obs=length(unique(obs$origin_gene))
    nbgenes.sim=length(unique(sim$origin_gene))

    print(paste(nrow(obs)," observed contacts for",nbgenes.obs,"genes after filtering target data"))
    print(paste(nrow(sim)," simulated contacts for",nbgenes.sim,"genes after filtering target data"))
        
    ## here we take only well-conserved enhancers!

    obs$align_score=pcungapped[obs$origin_enh]
    sim$align_score=pcungapped[sim$origin_enh]
    
    obs=obs[which(obs$align_score>=minAlignScore),]
    sim=sim[which(sim$align_score>=minAlignScore),]

    ## take also contacts that are in conserved synteny - gene and lifted enhancer on the same chromosome

    ## read synteny conservation
      
    synt_obs <- fread(paste(pathEvolution,"/synteny_conservation/", enh, "/", ref, "2", tg, "_original_synteny.txt", sep=""), header=T)
    synt_simul <- fread(paste(pathEvolution,"/synteny_conservation/", enh, "/", ref, "2", tg, "_simulated_synteny.txt", sep=""), header=T)
    
    class(synt_obs)<-"data.frame"
    class(synt_simul)<-"data.frame"
    
    synt_obs$id=paste(synt_obs$origin_gene, synt_obs$origin_enh, sep="-")
    synt_simul$id=paste(synt_simul$origin_gene, synt_simul$origin_enh, sep="-")

    synt_obs$target_dist=synt_obs$target_min_dist_allTSS
    synt_simul$target_dist=synt_simul$target_min_dist_allTSS

    ## if in trans, target distance is NA
    synt_obs$target_dist[which(synt_obs$target_dist=="trans")] <- NA
    synt_simul$target_dist[which(synt_simul$target_dist=="trans")] <- NA 
    
    ## transform distance to numeric values
    synt_obs$target_dist <- as.numeric(synt_obs$target_dist)
    synt_simul$target_dist <- as.numeric(synt_simul$target_dist)

    ## select interactions found between minDistance and maxDistance in tg species

    synt_obs=synt_obs[which(synt_obs$target_dist>=minDistance & synt_obs$target_dist<=maxDistance),]
    synt_simul=synt_simul[which(synt_simul$target_dist>=minDistance & synt_simul$target_dist<=maxDistance),]

    ## select interactions with conserved synteny

    obs=obs[which(obs$id%in%synt_obs$id),]
    sim=sim[which(sim$id%in%synt_simul$id),]
    
    ## save final data
    
    contact.conservation[[paste(ref, "2", tg, sep="")]][[enh]]=list("obs"=obs, "sim"=sim) 
    
  }
}

#######################################################################################

save(list=c("unfiltered.contact.conservation", "contact.conservation"), file=paste(pathFigures, "RData/data.contact.conservation.enhancers.RData",sep=""))

#######################################################################################

