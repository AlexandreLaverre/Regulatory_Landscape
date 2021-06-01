##########################################################################

library(Hmisc)

options(stringsAsFactors = FALSE)

#########################################################################

## if it's the first time we run this figure, we load and prepare data

objects=ls()

if(!"pathScripts"%in%objects){
  load=T
  prepare=T
  source("../main_figures/parameters.R")

  library(bootBCa, lib=pathRlibs)
  
  set.seed(19)
}

##################################################################

load(paste(pathFigures, "RData/data.sample.info.RData", sep=""))
load(paste(pathFigures, "RData/data.enhancer.statistics.RData", sep=""))
load(paste(pathFigures, "RData/data.contact.conservation.enhancers.RData", sep=""))

## prepare data for each species

chicago.dist <- list()
chicago.dist.conf.low <- list()
chicago.dist.conf.high <- list()

cons.seq <- list()
cons.seq.conf.low <- list()
cons.seq.conf.high <- list()

cons.contact <- list()
cons.contact.conf.low <- list()
cons.contact.conf.high <- list()

cons.synteny <- list()
cons.synteny.conf.low <- list()
cons.synteny.conf.high <- list()

nb_chicago_class = 10

for(sp in c("human", "mouse")){
  ## prepare data
  
  tg=setdiff(c("human", "mouse"), sp)
  sampleinfo.tg=sampleinfo[[tg]]
  
  chicago.dist[[sp]] <- list()
  chicago.dist.conf.low[[sp]] <- list()
  chicago.dist.conf.high[[sp]] <- list()
  
  cons.seq[[sp]] <- list()
  cons.seq.conf.low[[sp]] <- list()
  cons.seq.conf.high[[sp]] <- list()
  
  cons.contact[[sp]] <- list()
  cons.contact.conf.low[[sp]] <- list()
  cons.contact.conf.high[[sp]] <- list()
  
  cons.synteny[[sp]] <- list()
  cons.synteny.conf.low[[sp]] <- list()
  cons.synteny.conf.high[[sp]] <- list()
  
  load(paste(pathFigures, "RData/data.sequence.conservation.stats.pcungapped.", sp, ".RData", sep=""))
  load(paste(pathFigures, "RData/data.synteny.conservation.", sp,".RData", sep=""))
  
  for (enh in enhancer.datasets[[sp]]){
    ## Calculate sequence conservation 
    enh_obs <- enhancer.statistics[[sp]][[enh]][["original"]]
    enh_obs$class_score <- cut2(enh_obs$median_score, g=nb_chicago_class)
    
    enh_align <- list_align_enh[[enh]][["enh_align_obs"]]
    rownames(enh_align) <- enh_align$ID
    common=intersect(rownames(enh_obs), rownames(enh_align))
    
    enh_align = enh_align[common,]
    enh_align$class_score = enh_obs[common,]$class_score

    print(paste("confidence intervals, class score", enh, sp, tg))
    
    BC.seq=tapply(100*enh_align[[tg]], as.factor(enh_align$class_score), function(x) BCa(x, delta=NA, M=100, theta=mean, na.rm=T))
    cons.seq[[sp]][[enh]] <- unlist(lapply(BC.seq, function(x) x[3]))
    cons.seq.conf.low[[sp]][[enh]] <- unlist(lapply(BC.seq, function(x) x[4]))
    cons.seq.conf.high[[sp]][[enh]] <-unlist(lapply(BC.seq, function(x) x[5]))
    
    ## Calculate relation with distance to promoter
    contact_obs <- contact.conservation[[paste0(sp, "2", tg)]][[enh]][["obs"]] 
    
    contact_obs$class_dist <- cut(contact_obs$origin_dist, breaks=seq(from=minDistance, to=maxDistance, by=50000), include.lowest = T)
    
    print(paste("confidence intervals, chicago score", enh, sp, tg))
    
    BC.chic=tapply(contact_obs$origin_med_score, as.factor(contact_obs$class_dist), function(x) BCa(x, delta=NA, M=100, theta=mean, na.rm=T))
    chicago.dist[[sp]][[enh]] <- unlist(lapply(BC.chic, function(x) x[3]))
    chicago.dist.conf.low[[sp]][[enh]] <- unlist(lapply(BC.chic, function(x) x[4]))
    chicago.dist.conf.high[[sp]][[enh]] <-unlist(lapply(BC.chic, function(x) x[5]))
    
    ## Calculate contact conservation proportion
    contact_obs$cons=apply(contact_obs[,sampleinfo.tg$Sample.ID], 1, function(x) any(x>0))
    
    contact_obs$class_score <- cut2(contact_obs$origin_med_score, g=nb_chicago_class)
    
    cons.contact[[sp]][[enh]] <- sapply(levels(contact_obs$class_score), function(x)
                                        (nrow(contact_obs[which(contact_obs$cons == T & contact_obs$class_score == x),])/nrow(contact_obs[which(contact_obs$class_score == x ),]))*100)
    
    cons.contact.conf.low[[sp]][[enh]] <- sapply(levels(contact_obs$class_score), function(x)
                                                 prop.test(nrow(contact_obs[which(contact_obs$cons == T & contact_obs$class_score == x),]), nrow(contact_obs[which(contact_obs$class_score == x ),]))$conf.int[1]*100)
    
    cons.contact.conf.high[[sp]][[enh]] <- sapply(levels(contact_obs$class_score), function(x)
                                                  prop.test(nrow(contact_obs[which(contact_obs$cons == T & contact_obs$class_score == x),]), nrow(contact_obs[which(contact_obs$class_score == x ),]))$conf.int[2]*100)
    
    ## Calculate synteny conservation 
    synt_obs <- conserv_synteny[[enh]][[tg]][["synt_obs"]]
    rownames(synt_obs) <- paste0(synt_obs$origin_gene, "-", synt_obs$origin_enh)
    rownames(contact_obs) <- paste0(contact_obs$origin_gene, "-", contact_obs$origin_enh)
    
    common=intersect(rownames(synt_obs), rownames(contact_obs))
    synt_obs <- synt_obs[common,]
    synt_obs$class_score = contact_obs[common,]$class_score
    
    cons.synteny[[sp]][[enh]] <- sapply(levels(synt_obs$class_score), function(x)
                                        (nrow(synt_obs[which(synt_obs$target_dist <= maxDistanceSyntenyTarget & synt_obs$class_score == x),])/nrow(synt_obs[which(synt_obs$class_score == x ),]))*100)
    
    cons.synteny.conf.low[[sp]][[enh]] <- sapply(levels(contact_obs$class_score), function(x)
                                                 prop.test(nrow(synt_obs[which(synt_obs$target_dist <= maxDistanceSyntenyTarget & synt_obs$class_score == x),]), nrow(synt_obs[which(synt_obs$class_score == x ),]))$conf.int[1]*100)
    
    cons.synteny.conf.high[[sp]][[enh]] <- sapply(levels(contact_obs$class_score), function(x)
                                                  prop.test(nrow(synt_obs[which(synt_obs$target_dist <= maxDistanceSyntenyTarget & synt_obs$class_score == x),]), nrow(synt_obs[which(synt_obs$class_score == x ),]))$conf.int[2]*100)
    
  }
}

#####################################################################################

save(list=ls(), file=paste(pathFigures, "RData/data.bootstrap.chicago.scores.RData",sep=""))

#####################################################################################
