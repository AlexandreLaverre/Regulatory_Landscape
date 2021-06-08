####################################################################################################

library(data.table)

options(stringsAsFactors = FALSE)

source("parameters.R")

all.species <- c("human", "macaque", "dog", "cow", "elephant", "rabbit", "mouse", "rat", "opossum", "chicken")

load(paste(pathFigures, "RData/data.gene.enhancer.contacts.RData", sep=""))

####################################################################################################

for(ref_sp in c("human", "mouse")){
 
  target_sp = setdiff(c("human", "mouse"), ref_sp)

  species <- setdiff(all.species, ref_sp)
  
  ## define paths

  path_evol <- paste(pathFinalData, "SupplementaryDataset7/", ref_sp, "/", sep="")
  path_annot <- paste(pathFinalData, "SupplementaryDataset4/", ref_sp, "/", sep="")
  
  ## enhancer datasets

  enhancers <- enhancer.datasets[[ref_sp]]

  ## synteny conservation, for allspecies pairs
  
  conserv_synteny <- list()

  for (enh in enhancers){
    message("Running conserv synteny with ", enh)
    
    ## filtered gene enhancer contacts

    filtered.contacts.obs=gene.enhancer.contacts[[ref_sp]][[enh]][["real"]]
    filtered.contacts.sim=gene.enhancer.contacts[[ref_sp]][[enh]][["simulated"]]

    filtered.contacts.obs$id=paste(filtered.contacts.obs$gene, filtered.contacts.obs$enhancer, sep="-")
    filtered.contacts.sim$id=paste(filtered.contacts.sim$gene, filtered.contacts.sim$enhancer, sep="-")

    rownames(filtered.contacts.obs)=filtered.contacts.obs$id
    rownames(filtered.contacts.sim)=filtered.contacts.sim$id
    
    conserv_synteny[[enh]] <- list()
    
    for (sp in species){

      ## prepare result 
      conserv_synteny[[enh]][[sp]]<-list()
      
      ## load sequence conservation 
      load(paste(pathFigures, "RData/data.sequence.conservation.enhancers.",enh,".",ref_sp,"2", sp,".RData", sep=""))

      ## read synteny conservation
      
      synt_obs <- fread(paste(path_evol,"/synteny_conservation/", enh, "/", ref_sp, "2", sp, "_original_synteny.txt", sep=""), header=T)
      synt_simul <- fread(paste(path_evol,"/synteny_conservation/", enh, "/", ref_sp, "2", sp, "_simulated_synteny.txt", sep=""), header=T)

      class(synt_obs)<-"data.frame"
      class(synt_simul)<-"data.frame"

      print(paste(nrow(synt_obs)," observed contacts before filtering"))
      print(paste(nrow(synt_simul)," simulated contacts before filtering"))
      
      ## select only previously filtered gene - enhancer pairs

      synt_obs$id=paste(synt_obs$origin_gene, synt_obs$origin_enh, sep="-")
      synt_simul$id=paste(synt_simul$origin_gene, synt_simul$origin_enh, sep="-")

      synt_obs=synt_obs[which(synt_obs$id%in%filtered.contacts.obs$id),]
      synt_simul=synt_simul[which(synt_simul$id%in%filtered.contacts.sim$id),]

      print(paste(nrow(synt_obs)," observed contacts after filtering"))
      print(paste(nrow(synt_simul)," simulated contacts after filtering"))

     
               
      ## threshold alignment score: 5% quantile, observed values
            
      synt_obs$align_score=pcungapped[synt_obs$origin_enh]
      synt_simul$align_score=pcungapped[synt_simul$origin_enh]
    
      align.threshold=quantile(synt_obs$align_score, p=0.1, na.rm=T)

      print(paste("alignment score threshold", align.threshold))
      
      ## enhancers aligned at least at align.threshold in target genome
      synt_obs <- synt_obs[which(synt_obs$align_score >= align.threshold),]
      synt_simul <- synt_simul[which(synt_simul$align_score >=align.threshold),]

      ## we use the minimum distance between baited TSS and enhancer for the reference species
      ## the minimum distance between all TSS and enhancer for target species

      synt_obs$origin_dist=synt_obs$origin_min_dist_baitedTSS
      synt_simul$origin_dist=synt_simul$origin_min_dist_baitedTSS

      synt_obs$target_dist=synt_obs$target_min_dist_allTSS
      synt_simul$target_dist=synt_simul$target_min_dist_allTSS

      ## compute distance classes for original pairs

      synt_obs$class_dist <- cut(synt_obs$origin_dist, breaks=seq(from=minDistance, to=maxDistance, by=50000), include.lowest = T)
      synt_simul$class_dist <- cut(synt_simul$origin_dist, breaks=seq(from=minDistance, to=maxDistance, by=50000), include.lowest = T)
      
      ## if in trans, target distance is NA
      synt_obs$target_dist[which(synt_obs$target_dist=="trans")] <- NA
      synt_simul$target_dist[which(synt_simul$target_dist=="trans")] <- NA 

      ## transform distance to numeric values
      synt_obs$target_dist <- as.numeric(synt_obs$target_dist)
      synt_simul$target_dist <- as.numeric(synt_simul$target_dist)

      ## remove gene coordinates etc

      synt_obs=synt_obs[,which(!(colnames(synt_obs)%in%c("origin_gene_coord", "target_gene_coord", "nb_genes_inbetween", "origin_min_dist_baitedTSS", "origin_min_dist_allTSS", "target_min_dist_allTSS")))]

      synt_simul=synt_simul[,which(!(colnames(synt_simul)%in%c("origin_gene_coord", "target_gene_coord", "nb_genes_inbetween", "origin_min_dist_baitedTSS", "origin_min_dist_allTSS", "target_min_dist_allTSS")))]

      ## add median score for contacts

      synt_obs$median_score=filtered.contacts.obs[synt_obs$id, "median_score"]
      synt_simul$median_score=filtered.contacts.sim[synt_simul$id, "median_score"]
      
      
      conserv_synteny[[enh]][[sp]]=list("synt_obs"=synt_obs, "synt_simul"=synt_simul)
            
      message(sp, " done !")
    }
  }

  ## save data 
  
  save(conserv_synteny, file = paste(pathFigures, "/RData/data.synteny.conservation.", ref_sp, ".RData", sep=""))
  
  ## all done!
}

#######################################################################################################################################
