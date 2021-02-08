##############################################################################

source("parameters.R")
library(Hmisc)

##############################################################################

load(paste(pathFigures, "RData/data.sample.info.RData", sep=""))
load(paste(pathFigures, "RData/data.contact.conservation.enhancers.RData", sep=""))

##############################################################################

gene.ontologies = c("all", "dvpt", "immune", "other")

##############################################################################

for(ref in c("human", "mouse")){
  
  tg=setdiff(c("human", "mouse"), ref)
  
  enhancers=enhancer.datasets[[ref]]
  contact.cons=contact.conservation[[paste0(ref, "2", tg)]] ## contact conservation for enhancers, already filtered
  
  sampleinfo.ref=sampleinfo[[ref]]
  sampleinfo.tg=sampleinfo[[tg]]

  ## select gene types
  
  dvpt = unlist(read.table(paste(pathFinalData, "SupplementaryDataset3/gene_ontology/", ref, "_developmental_process_and_child_terms_genes_Ensembl94.txt", sep="")), use.names = FALSE)
  immune = unlist(read.table(paste(pathFinalData, "SupplementaryDataset3/gene_ontology/", ref, "_immune_process_and_child_terms_genes_Ensembl94.txt", sep="")), use.names = FALSE)
  
  genes <- list("dvpt"=dvpt, "immune"=immune)
  
  ##############################################################################
  
  cons.stats <- list()
  cons.dist <- list()
  cons.dist.conf.low <- list()
  cons.dist.conf.high <- list()
  cons.nb.cell <- list()
  cons.nb.cell.conf.low <- list()
  cons.nb.cell.conf.high <- list()
  cons.common.cell <- list()
  cons.common.cell.conf.low <- list()
  cons.common.cell.conf.high <- list()
    
  for(enh in enhancers){

    cc.obs <- contact.cons[[enh]][["obsobs"]]
    cc.sim <- contact.cons[[enh]][["simsim"]]
    
    cc.obs$cons <- apply(cc.obs[,sampleinfo.tg$Sample.ID], 1, function(x) any(x>0))
    cc.sim$cons <- apply(cc.sim[,sampleinfo.tg$Sample.ID], 1, function(x) any(x>0))
    
    for (onto in gene.ontologies){
      if (onto == "all"){
        obs <- cc.obs
        sim <- cc.sim
      } else if(onto == "other"){
        obs <- cc.obs[which(cc.obs$origin_gene %nin% genes[["dvpt"]] & cc.obs$origin_gene %nin% genes[["immune"]]),]
        sim <- cc.sim[which(cc.sim$origin_gene %nin% genes[["dvpt"]] & cc.sim$origin_gene %nin% genes[["immune"]]),]
      }else{
        obs <- cc.obs[which(cc.obs$origin_gene %in% genes[[onto]]),]
        sim <- cc.sim[which(cc.sim$origin_gene %in% genes[[onto]]),]
      }
      
      pc.cons.obs <- 100*length(which(obs$cons))/dim(obs)[1]
      pc.cons.sim <- 100*length(which(sim$cons))/dim(sim)[1]
      
      test.obs <- prop.test(length(which(obs$cons)), dim(obs)[1])
      test.sim <- prop.test(length(which(sim$cons)), dim(sim)[1])
      
      cons.stats[[enh]][[onto]]=list("pc.cons.obs"=pc.cons.obs, "pc.cons.sim"=pc.cons.sim,
                         "test.obs"=test.obs, "test.sim"=test.sim)
    }
    
     ################ contact conservation in common cell types   ################ 
    
    cells <- c("ESC", "preadip", "Bcell")

    common.cell<-list()
    
    for(c in cells){
      common.cell[[c]]=list()
    }
    
    for(sp in c(ref, tg)){
      this.sampleinfo<-sampleinfo[[sp]]
      common.cell[["ESC"]][[sp]]<-this.sampleinfo$Sample.ID[which(this.sampleinfo[,"Broad.cell.type.or.tissue"]=="embryonic stem cells")]
      common.cell[["preadip"]][[sp]]<-this.sampleinfo$Sample.ID[which(this.sampleinfo[,"Broad.cell.type.or.tissue"]=="pre-adipocytes")]
      common.cell[["Bcell"]][[sp]]<-this.sampleinfo$Sample.ID[which(this.sampleinfo[,"Broad.cell.type.or.tissue"]=="B lymphocytes")]
    }
    
    cons.stats.cell <- list()
    
    for (cell in cells){
      ## select enhancers that are contacted in this cell type
      
      cc.obs.cell <- cc.obs[which(apply(cc.obs[common.cell[[cell]][[ref]]], 1, function(x) any(x>0))),]
      cc.sim.cell <- cc.sim[which(apply(cc.sim[common.cell[[cell]][[ref]]], 1, function(x) any(x>0))),]
      
      cc.obs.cell$cons<-apply(cc.obs.cell[common.cell[[cell]][[tg]]], 1, function(x) any(x>0))
      cc.sim.cell$cons<-apply(cc.sim.cell[common.cell[[cell]][[tg]]], 1, function(x) any(x>0))
      
      pc.cons.obs <- 100*length(which(cc.obs.cell$cons))/dim(cc.obs.cell)[1]
      pc.cons.sim <- 100*length(which(cc.sim.cell$cons))/dim(cc.sim.cell)[1]
      
      test.obs <- prop.test(length(which(cc.obs.cell$cons)), dim(cc.obs.cell)[1])
      test.sim <- prop.test(length(which(cc.sim.cell$cons)), dim(cc.sim.cell)[1])
      
      cons.stats.cell[[cell]]=list("pc.cons.obs" = pc.cons.obs, "pc.cons.sim" = pc.cons.sim, "test.obs" = test.obs, "test.sim" = test.sim)
    }
    
    cons.common.cell[[enh]] = sapply(cells, function(x) c(cons.stats.cell[[x]]$pc.cons.obs, cons.stats.cell[[x]]$pc.cons.sim))
    cons.common.cell.conf.low[[enh]] =  sapply(cells, function(x) c(cons.stats.cell[[x]]$test.obs$conf.int[1]*100, cons.stats.cell[[x]]$test.sim$conf.int[1]*100))
    cons.common.cell.conf.high[[enh]] =  sapply(cells, function(x) c(cons.stats.cell[[x]]$test.obs$conf.int[2]*100, cons.stats.cell[[x]]$test.sim$conf.int[2]*100))

   ############### contact conservation according to distance class   ################ 
   ## minDistance and maxDistance defined in parameters.R
    
    cc.obs$dist_class <- cut(cc.obs$origin_dist, breaks=seq(from=minDistance, to=maxDistance, by=50e3), include.lowest = T)
    cc.sim$dist_class <- cut(cc.sim$origin_dist, breaks=seq(from=minDistance, to=maxDistance, by=50e3), include.lowest = T)
    
    for (onto in gene.ontologies){
      if (onto == "all"){
        obs = cc.obs
        sim = cc.sim
      }else if(onto == "other"){
        obs = cc.obs[which(cc.obs$origin_gene %nin% genes[["dvpt"]] & cc.obs$origin_gene %nin% genes[["immune"]]),]
        sim = cc.sim[which(cc.sim$origin_gene %nin% genes[["dvpt"]] & cc.sim$origin_gene %nin% genes[["immune"]]),]
        
      }else{
        obs = cc.obs[which(cc.obs$origin_gene %in% genes[[onto]]),]
        sim = cc.sim[which(cc.sim$origin_gene %in% genes[[onto]]),]
      }
      
      contact.data <-  list("obs"=obs, "sim"=sim)
      
      cons.dist[[enh]][[onto]] <- t(sapply(contact.data, function(y) sapply(levels(y$dist_class), function(x)
                                                                            (nrow(y[which(y$cons == T & y$dist_class == x),])/nrow(y[which(y$dist_class == x ),]))*100)))
      
      cons.dist.conf.low[[enh]][[onto]] <- t(sapply(contact.data, function(y) sapply(levels(y$dist_class), function(x)
        prop.test(nrow(y[which(y$cons == T & y$dist_class == x),]), nrow(y[which(y$dist_class == x ),]))$conf.int[1]*100)))
      
      cons.dist.conf.high[[enh]][[onto]] <- t(sapply(contact.data, function(y) sapply(levels(y$dist_class), function(x)
        prop.test(nrow(y[which(y$cons == T & y$dist_class == x),]), nrow(y[which(y$dist_class == x ),]))$conf.int[2]*100)))
      
    }    

    ################ contact conservation according to number of cell types   ################
    
    samples=sampleinfo.ref$Sample.ID 
    celltypes=sampleinfo.ref$Broad.cell.type.or.tissue
    names(celltypes)=samples
    
    if (ref == "human"){
      max.nb.cell = 7
    } else{
      max.nb.cell = 5
    }
    
    cc.obs$nb_celltypes <- apply(cc.obs[,samples],1, function(x) length(unique(celltypes[which(x > 0)])))
    cc.obs$celltype_class<- cut(cc.obs$nb_celltypes, breaks=c(0:max.nb.cell, max(cc.obs$nb_celltypes)), include.lowest=T)
    levels(cc.obs$celltype_class)=c(as.character(1:max.nb.cell), paste0(">", max.nb.cell))
    
    cc.sim$nb_celltypes <- apply(cc.sim[,samples],1, function(x) length(unique(celltypes[which(x > 0)])))
    cc.sim$celltype_class<- cut(cc.sim$nb_celltypes, breaks=c(0:max.nb.cell, max(cc.sim$nb_celltypes)), include.lowest=T)
    levels(cc.sim$celltype_class)=c(as.character(1:max.nb.cell), paste0(">", max.nb.cell))
    
    contact.data <-  list("obs"=cc.obs, "sim"=cc.sim)
    
    cons.nb.cell[[enh]] <- t(sapply(contact.data, function(y) sapply(levels(y$celltype_class), function(x)
      (nrow(y[which(y$cons == T & y$celltype_class == x),])/nrow(y[which(y$celltype_class == x ),]))*100)))
    
    cons.nb.cell.conf.low[[enh]] <- t(sapply(contact.data, function(y) sapply(levels(y$celltype_class), function(x)
      prop.test(nrow(y[which(y$cons == T & y$celltype_class == x),]), nrow(y[which(y$celltype_class == x ),]))$conf.int[1]*100)))
    
    cons.nb.cell.conf.high[[enh]] <- t(sapply(contact.data, function(y) sapply(levels(y$celltype_class), function(x)
      prop.test(nrow(y[which(y$cons == T & y$celltype_class == x),]), nrow(y[which(y$celltype_class == x ),]))$conf.int[2]*100)))
    
  }
  
  cons <- list()
  cons.conf.low <- list()
  cons.conf.high <- list()
  
  for (onto in gene.ontologies){
    cons[[onto]] = sapply(enhancers, function(x) c(cons.stats[[x]][[onto]]$pc.cons.obs, cons.stats[[x]][[onto]]$pc.cons.sim))
    cons.conf.low[[onto]] =  sapply(enhancers, function(x) c(cons.stats[[x]][[onto]]$test.obs$conf.int[1]*100, cons.stats[[x]][[onto]]$test.sim$conf.int[1]*100))
    cons.conf.high[[onto]] =  sapply(enhancers, function(x) c(cons.stats[[x]][[onto]]$test.obs$conf.int[2]*100, cons.stats[[x]][[onto]]$test.sim$conf.int[2]*100))
  }
  
  
###################### output ######################
  
  save(cons, cons.conf.low, cons.conf.high,
       cons.dist, cons.dist.conf.low, cons.dist.conf.high,
       cons.nb.cell, cons.nb.cell.conf.low, cons.nb.cell.conf.high,
       cons.common.cell, cons.common.cell.conf.low, cons.common.cell.conf.high, 
       file = paste(pathFigures, "RData/data.contact.conservation.enhancers.", ref, ".stats.RData", sep=""))
  
}


#######################################################################################
