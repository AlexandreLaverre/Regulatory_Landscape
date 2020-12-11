#######################################################################################
library(data.table)
library(Hmisc)
options(stringsAsFactors = FALSE)

source("parameters.R")

pathEvolution=paste(pathFinalData, "SupplementaryDataset7", sep="")


#######################################################################################

load(paste(pathFigures,"RData/data.enhancer.statistics.RData", sep=""))
load(paste(pathFigures, "RData/data.sample.info.RData", sep=""))

#######################################################################################

genes.conservation=list()

for(ref in c("human", "mouse")){
  
  tg=setdiff(c("human", "mouse"), ref)

  sampleinfo.ref=sampleinfo[[ref]]
  sampleinfo.tg=sampleinfo[[tg]]
  
  cells <- c("ESC", "preadip", "Bcell")
  
  ESC.common = list("human"= c("hESC"), "mouse"=c("ESC", "ESC_18", "ESC_wild"))
  preadip.common = list("human"= c("pre_adipo"), "mouse"=c("preadip_D0", "preadip_D2", "preadip_4H"))
  Bcell.common = list("human"= c("TB", "NB"), "mouse"=c("preB_aged", "preB_young"))
  common.cell <- list("ESC"=ESC.common, "preadip"=preadip.common, "Bcell"=Bcell.common)
  
  
  for(enh in enhancer.datasets[[ref]]){
    print(enh)
    
    obsobs=fread(paste(pathEvolution, "/", ref, "/contact_conservation/", enh, "/", ref, "_original2", tg,"_original.txt_test", sep=""), h=T, stringsAsFactors=F, sep="\t")
    simsim=fread(paste(pathEvolution, "/", ref, "/contact_conservation/", enh, "/", ref, "_simulated2", tg,"_simulated.txt_test", sep=""), h=T, stringsAsFactors=F, sep="\t")
    
    class(obsobs)="data.frame"
    class(simsim)="data.frame"
    
    ## gene with orthologue one2one in target specie
    obsobs = obsobs[which(!is.na(obsobs$target_gene)),]
    simsim = simsim[which(!is.na(simsim$target_gene)),]
    
    ## filtered enhancers
    obsobs=obsobs[which(obsobs$BLAT_match == 1 & obsobs$origin_dist >= minDistance),]
    simsim=simsim[which(simsim$BLAT_match == 1 & simsim$origin_dist >= minDistance),]
    
    # Class of distances
    obsobs$class_dist = cut(obsobs$origin_dist, breaks=c(minDistance, 100000, 500000, maxDistance), include.lowest = T)
    simsim$class_dist = cut(simsim$origin_dist, breaks=c(minDistance, 100000, 500000, maxDistance), include.lowest = T)

    # Conservation of contact
    obsobs$cons=apply(obsobs[,sampleinfo.tg$Sample.ID], 1, function(x) any(x>0))
    simsim$cons=apply(simsim[,sampleinfo.tg$Sample.ID], 1, function(x) any(x>0))
    
    datas <- list("obs"=obsobs, "sim"=simsim)
    
    for (data.name in c("obs", "sim")){
      print(data.name)
      
      data = datas[[data.name]]
      all_genes = levels(factor(data$origin_gene))
  
      # Calculate conservation for each distance class
      for (dist in c(levels(data$class_dist), "all")){
        
        print(dist)
        if (dist == "all"){selected_dist = data}else{selected_dist = data[which(data$class_dist == dist),]}

        nb_total = unlist(with(selected_dist, tapply(origin_enh, factor(origin_gene, levels=all_genes), function(x) length(x))))
        align_score =  unlist(with(selected_dist, tapply(align_score, factor(origin_gene, levels=all_genes), median, na.rm=T)))
        seq_conserv = unlist(with(selected_dist, tapply(align_score, factor(origin_gene, levels=all_genes), function(x) length(which(x>=0.4)))))
        
        selected_align = selected_dist[which(selected_dist$align_score >= 0.4),]
        synt_conserv = with(selected_align, tapply(target_dist, factor(origin_gene, levels=all_genes), function(x) length(which(x <= maxDistanceSyntenyTarget))))
        
        selected_gene = selected_align[which(selected_align$target_data == TRUE),] # ortologous genes present in PCHIC in target specie
        contact_conserv = with(selected_gene, tapply(cons, factor(origin_gene, levels=all_genes), function(x) length(which(x == TRUE))))
      
        genes.conservation[[enh]][[data.name]][[dist]] <- data.frame("nb_total"=nb_total, "align_score"=align_score, "seq_conserv"=seq_conserv, 
                                                                     "synt_conserv"=synt_conserv, "contact_conserv"=contact_conserv)
        # Calculate ratio with constraints on minimum number
        if (enh == "FANTOM5"){nb_min = 2}else{nb_min=5}
        
        genes.conservation[[enh]][[data.name]][[dist]]$seq_conserv <- with(genes.conservation[[enh]][[data.name]][[dist]], ifelse(nb_total >= nb_min | nb_total <= 100,  seq_conserv, NA))
        genes.conservation[[enh]][[data.name]][[dist]]$ratio_cons_seq <- with(genes.conservation[[enh]][[data.name]][[dist]], ifelse(nb_total >= nb_min & nb_total <= 100, seq_conserv/nb_total, NA))
        genes.conservation[[enh]][[data.name]][[dist]]$ratio_cons_synt <- with(genes.conservation[[enh]][[data.name]][[dist]], ifelse(seq_conserv >= nb_min & seq_conserv <= 100,  synt_conserv/seq_conserv, NA))
        genes.conservation[[enh]][[data.name]][[dist]]$ratio_cons_int <- with(genes.conservation[[enh]][[data.name]][[dist]], ifelse(seq_conserv >= nb_min & seq_conserv <= 100,  contact_conserv/seq_conserv, NA))
        
        # Class of conservation ratio
        genes.conservation[[enh]][[data.name]][[dist]]$class_nb_contact = cut2(genes.conservation[[enh]][[data.name]][[dist]]$nb_total, g=5, include.lowest=T)
        genes.conservation[[enh]][[data.name]][[dist]]$class_cons_seq = cut(genes.conservation[[enh]][[data.name]][[dist]]$ratio_cons_seq, breaks=c(0, 0.10, 0.25, 0.5, 0.75, 1), include.lowest=T)
        genes.conservation[[enh]][[data.name]][[dist]]$class_cons_synt = cut(genes.conservation[[enh]][[data.name]][[dist]]$ratio_cons_synt, breaks=c(0, 0.99, 1), include.lowest=T)
        genes.conservation[[enh]][[data.name]][[dist]]$class_cons_cont = cut(genes.conservation[[enh]][[data.name]][[dist]]$ratio_cons_int,  breaks=c(0, 0.01, 0.25, 0.5, 0.75, 1), include.lowest=T)
        genes.conservation[[enh]][[data.name]][[dist]]$class_align_score = cut2(genes.conservation[[enh]][[data.name]][[dist]]$align_score, g=5, include.lowest=T)
        
      }
      names(genes.conservation[[enh]][[data.name]]) = c("25kb - 100kb", "100kb - 500kb", "500kb - 2Mb", "all")
    }
  }
  
  ## save data
  save(genes.conservation, file=paste(pathFigures, "RData/data.", ref, ".gene.regland.conservation.RData",sep=""))
}


