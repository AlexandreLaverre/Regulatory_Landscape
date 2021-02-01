#######################################################################################
library(data.table)
library(Hmisc)
options(stringsAsFactors = FALSE)

source("parameters.R")

pathEvolution=paste(pathFinalData, "SupplementaryDataset7", sep="")

#######################################################################################

load(paste(pathFigures,"RData/data.enhancer.statistics.RData", sep=""))
load(paste(pathFigures, "RData/data.contact.conservation.enhancers.RData",sep=""))
load(paste(pathFigures, "RData/data.sample.info.RData", sep=""))
load(paste(pathFigures, "RData/data.ortho.genes.RData", sep=""))

#######################################################################################

genes.conservation=list()

for(ref in c("human", "mouse")){
  
  tg=setdiff(c("human", "mouse"), ref)

  sampleinfo.ref=sampleinfo[[ref]]
  sampleinfo.tg=sampleinfo[[tg]]
  
  for(enh in enhancer.datasets[[ref]]){
    print(enh)
    
    obs=contact.conservation[[paste(ref, "2", tg, sep="")]][[enh]][["obs"]]
    sim=contact.conservation[[paste(ref, "2", tg, sep="")]][[enh]][["sim"]]

    ## select previously filtered ortho genes

    obs=obs[which(obs$origin_gene%in%ortho[,ref] & obs$target_gene%in%ortho[,tg]),]
    sim=sim[which(sim$origin_gene%in%ortho[,ref] & sim$target_gene%in%ortho[,tg]),]

    ## gene must have bait in other species

     obs=obs[which(obs$target_data == "TRUE"),]
    sim=sim[which(sim$target_data == "TRUE"),]
        
    
    ## Class of distances
    obs$class_dist = cut(obs$origin_dist, breaks=c(minDistance, 100000, 500000, maxDistance), include.lowest = T)
    sim$class_dist = cut(sim$origin_dist, breaks=c(minDistance, 100000, 500000, maxDistance), include.lowest = T)

    ## Conservation of contact
    obs$cons=apply(obs[,sampleinfo.tg$Sample.ID], 1, function(x) any(x>0))
    sim$cons=apply(sim[,sampleinfo.tg$Sample.ID], 1, function(x) any(x>0))
    
    data.list <- list("obs"=obs, "sim"=sim)
    
    for (data.name in c("obs", "sim")){
      print(data.name)
      
      data = data.list[[data.name]]
      all_genes = unique(data$origin_gene)
  
      # Calculate conservation for each distance class
      for (dist in c(levels(data$class_dist), "all")){
        
        print(dist)
        if (dist == "all"){
          selected_dist = data
        } else{
          selected_dist = data[which(data$class_dist == dist),]
        }

        nb_total = unlist(with(selected_dist, tapply(origin_enh, factor(origin_gene, levels=all_genes), function(x) length(x))))
        align_score =  unlist(with(selected_dist, tapply(align_score, factor(origin_gene, levels=all_genes), median, na.rm=T)))
        seq_conserv = unlist(with(selected_dist, tapply(align_score, factor(origin_gene, levels=all_genes), function(x) length(which(x>= minAlignScore)))))
        
        selected_align = selected_dist[which(selected_dist$align_score >= minAlignScore),]
        synt_conserv = with(selected_align, tapply(target_dist, factor(origin_gene, levels=all_genes), function(x) length(which(as.numeric(x) <= maxDistanceSyntenyTarget))))
        
        selected_gene = selected_align[which(selected_align$target_data == TRUE),] # ortologous genes present in PCHIC in target specie
        contact_conserv = with(selected_gene, tapply(cons, factor(origin_gene, levels=all_genes), function(x) length(which(x == TRUE))))
      
        genes.conservation[[enh]][[data.name]][[dist]] <- data.frame("nb_total"=nb_total, "align_score"=align_score, "seq_conserv"=seq_conserv, 
                                                                     "synt_conserv"=synt_conserv, "contact_conserv"=contact_conserv)
        ## Calculate ratio with constraints on minimum number of contacted enhancers per gene
        
        if (enh == "FANTOM5"){
          nb_min = 2
        }else{
          nb_min=5
        }
        
        genes.conservation[[enh]][[data.name]][[dist]]$seq_conserv <- with(genes.conservation[[enh]][[data.name]][[dist]], ifelse(nb_total >= nb_min | nb_total <= 100,  seq_conserv, NA))
        genes.conservation[[enh]][[data.name]][[dist]]$ratio_cons_seq <- with(genes.conservation[[enh]][[data.name]][[dist]], ifelse(nb_total >= nb_min & nb_total <= 100, seq_conserv/nb_total, NA))
        genes.conservation[[enh]][[data.name]][[dist]]$ratio_cons_synt <- with(genes.conservation[[enh]][[data.name]][[dist]], ifelse(seq_conserv >= nb_min & seq_conserv <= 100,  synt_conserv/seq_conserv, NA))
        genes.conservation[[enh]][[data.name]][[dist]]$ratio_cons_int <- with(genes.conservation[[enh]][[data.name]][[dist]], ifelse(seq_conserv >= nb_min & seq_conserv <= 100,  contact_conserv/seq_conserv, NA))
        
        # Class of conservation ratio
        genes.conservation[[enh]][[data.name]][[dist]]$class_nb_contact = cut2(genes.conservation[[enh]][[data.name]][[dist]]$nb_total, g=5, include.lowest=T)
        genes.conservation[[enh]][[data.name]][[dist]]$class_cons_seq = cut(genes.conservation[[enh]][[data.name]][[dist]]$ratio_cons_seq, breaks=c(0, 0.10, 0.25, 0.5, 0.75, 1), include.lowest=T)
        genes.conservation[[enh]][[data.name]][[dist]]$class_cons_synt = cut(genes.conservation[[enh]][[data.name]][[dist]]$ratio_cons_synt, breaks=c(0, 0.99, 1), include.lowest=T)
        genes.conservation[[enh]][[data.name]][[dist]]$class_cons_cont = cut(genes.conservation[[enh]][[data.name]][[dist]]$ratio_cons_int,  breaks=c(0,  0.5, 1), include.lowest=T)
        genes.conservation[[enh]][[data.name]][[dist]]$class_align_score = cut2(genes.conservation[[enh]][[data.name]][[dist]]$align_score, g=5, include.lowest=T)
        
      }
      names(genes.conservation[[enh]][[data.name]]) = c("25kb - 100kb", "100kb - 500kb", "500kb - 2Mb", "all")
    }
  }
  
  ## save data
  save(genes.conservation, file=paste(pathFigures, "RData/data.", ref, ".regland.conservation.RData",sep=""))
}

###########################################################################
