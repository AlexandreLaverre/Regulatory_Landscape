################################################################################################################################################
library(data.table)
options(stringsAsFactors = FALSE)

#setwd("/home/laverre/Data/Regulatory_landscape/scripts/main_figures")
#########################################################################

## if it's the first time we run this figure, we load and prepare data

objects=ls()

if(!"pathScripts"%in%objects){
  load=T
  prepare=T
  source("../main_figures/parameters.R")
}

##################################################################

if(load){
  sp="human"
  
  # Choose genes within : all ; dvpt ; other
  selected_genes = "all" 
  if (selected_genes == "dvpt"){`%get%` <- `%in%`}else{`%get%` <- Negate(`%in%`)}

  minDistance=25e3
  maxDistance=2e6
  
  load(paste(pathFigures, "RData/data.gene.enhancer.contacts.RData", sep=""))
  load(paste(pathFigures, "RData/data.gene.enhancer.contacts.conservation.RData", sep=""))
  load(paste(pathFigures, "RData/data.enhancer.statistics.RData", sep=""))
  
  path_overlap=paste(pathFinalData, "SupplementaryDataset7/", sp, "/sequence_conservation/enhancers/", sep="")
  
  if(selected_genes != "all"){
    dev_gene <- read.table(paste(pathFinalData, "SupplementaryDataset3/genes/", sp, "_dvpt_process_genes.txt", sep=""), header=T, sep="\t")
  }
    
  load=FALSE
}


##############################################################################################################################
################################################# Conserved contact global   #################################################
data_global <- c()
conf_low_global <- c()
conf_up_global <- c()
n_total_global <- c()

for (enh in enhancer.datasets[[sp]]){
  message(enh)
  all_obs <- gene.enhancer.contacts[[sp]][[enh]][["real"]] 
  all_simul <- gene.enhancer.contacts[[sp]][[enh]][["simulated"]] 
  
  contact_obs <- contacts.conservation[[sp]][[enh]][["real"]] 
  contact_simul <- contacts.conservation[[sp]][[enh]][["simulated"]] 
  
  # Select specific genes
  if(selected_genes != "all"){
    all_obs <- all_obs[which(all_obs$gene %get% dev_gene$Gene.ID),]
    all_simul <- all_simul[which(all_simul$gene %get% dev_gene$Gene.ID),]
    contact_obs <- contact_obs[which(contact_obs$origin_gene %get% dev_gene$Gene.ID),]
    contact_simul <- contact_simul[which(contact_simul$origin_gene %get% dev_gene$Gene.ID),]
  }
  
  
  # Select only contacts with no duplicated enhancers
  obs_stats <- enhancer.statistics[[sp]][[enh]][["real"]]
  simul_stats <- enhancer.statistics[[sp]][[enh]][["simulated"]]
  obs_stats$enh <-  do.call(paste,c(obs_stats[c("chr","start","end")],sep=":"))
  simul_stats$enh <-  do.call(paste,c(simul_stats[c("chr","start","end")],sep=":"))
  obs_stats <- obs_stats[which(obs_stats$BLAT_match < 2),]
  simul_stats <- simul_stats[which(simul_stats$BLAT_match < 2),]

  all_obs <- all_obs[which(all_obs$enhancer %in% obs_stats$enh),]
  all_simul <- all_simul[which(all_simul$enhancer %in% simul_stats$enh),]
  contact_obs <- contact_obs[which(contact_obs$origin_enh %in% obs_stats$enh),]
  contact_simul <- contact_simul[which(contact_simul$origin_enh %in% simul_stats$enh),]
  
  # Calculate proportion
  mat <- matrix(c(nrow(contact_obs), nrow(contact_simul), nrow(all_obs)-nrow(contact_obs), nrow(all_simul)-nrow(contact_simul)),2)
  
  conf_low_global <- append(conf_low_global, c((prop.test(x = nrow(contact_obs), n=nrow(all_obs), p=0.5)$conf.int[1])*100,
                                 (prop.test(x = nrow(contact_simul), n=nrow(all_simul), p=0.5)$conf.int[1])*100))
  
  conf_up_global <- append(conf_up_global, c((prop.test(x = nrow(contact_obs), n=nrow(all_obs), p=0.5)$conf.int[2])*100,
                               (prop.test(x = nrow(contact_simul), n=nrow(all_simul), p=0.5)$conf.int[2])*100))
  
  data_global <- append(data_global, c( (nrow(contact_obs)/nrow(all_obs))*100, (nrow(contact_simul)/nrow(all_simul))*100))
  
  n_total_global <- append(n_total_global, c(paste0("N = ", nrow(all_obs)), paste0("N = ", nrow(all_simul))))
  
  # Overlap with target enhancer
  if (enh %in% c("FANTOM5", "ENCODE")){
    overlap_enh <- read.table(paste(path_overlap, enh, "/enhancer_overlap_target_enhancer.txt", sep=""), header=T, sep="\t")
    enh_tot = nrow(overlap_enh)
    overlap_enh <- overlap_enh[which(overlap_enh[,5] != "NA"),]
    enh_overlap = nrow(overlap_enh)
    message("Proportion of lifted_enh overlap target enh : ", enh_overlap, " on ", enh_tot, " = ", enh_overlap/enh_tot  )
    
    all_obs <- all_obs[which(all_obs$enhancer %in% overlap_enh[,1]),]
    all_simul <- all_simul[which(all_simul$enhancer %in% overlap_enh[,1]),]
    contact_obs <- contact_obs[which(contact_obs$origin_enh %in% overlap_enh[,1]),]
    contact_simul <- contact_simul[which(contact_simul$origin_enh %in% overlap_enh[,1]),]
    mat <- matrix(c(nrow(contact_obs), nrow(contact_simul), nrow(all_obs)-nrow(contact_obs), nrow(all_simul)-nrow(contact_simul)),2)
    
    conf_low_global <- append(conf_low_global, c((prop.test(x = nrow(contact_obs), n=nrow(all_obs), p=0.5)$conf.int[1])*100,
                                                 (prop.test(x = nrow(contact_simul), n=nrow(all_simul), p=0.5)$conf.int[1])*100))
    
    conf_up_global <- append(conf_up_global, c((prop.test(x = nrow(contact_obs), n=nrow(all_obs), p=0.5)$conf.int[2])*100,
                                               (prop.test(x = nrow(contact_simul), n=nrow(all_simul), p=0.5)$conf.int[2])*100))
    
    data_global <- append(data_global, c( (nrow(contact_obs)/nrow(all_obs))*100, (nrow(contact_simul)/nrow(all_simul))*100))
    
    n_total_global <- append(n_total_global, c(paste0("N = ", nrow(all_obs)), paste0("N = ", nrow(all_simul))))
    
  }
  
  conf_low_global <- append(conf_low_global, NA)
  conf_up_global <- append(conf_up_global, NA)
  data_global <- append(data_global, NA)
  n_total_global <- append(n_total_global, NA)
}

conserv <- data.frame(data=data_global, conf_low=conf_low_global, conf_up=conf_up_global, n_total=n_total_global)
assign("conserv_global", conserv)

################################################# Conserv ~ genomic distance   #################################################

for (enh in enhancer.datasets[[sp]]){
  all_obs <- gene.enhancer.contacts[[sp]][[enh]][["real"]] 
  all_simul <- gene.enhancer.contacts[[sp]][[enh]][["simulated"]] 
  
  contact_obs <- contacts.conservation[[sp]][[enh]][["real"]] 
  contact_simul <- contacts.conservation[[sp]][[enh]][["simulated"]] 
  
  # Select specific genes
  if(selected_genes != "all"){
    all_obs <- all_obs[which(all_obs$gene %get% dev_gene$Gene.ID),]
    all_simul <- all_simul[which(all_simul$gene %get% dev_gene$Gene.ID),]
    contact_obs <- contact_obs[which(contact_obs$origin_gene %get% dev_gene$Gene.ID),]
    contact_simul <- contact_simul[which(contact_simul$origin_gene %get% dev_gene$Gene.ID),]
  }

  # Distance classes
  all_obs$class_dist <-cut(all_obs$dist, breaks=seq(from=0, to=maxDistance, by=50000), include.lowest = T)
  all_simul$class_dist <-cut(all_simul$dist, breaks=seq(from=0, to=maxDistance, by=50000), include.lowest = T)
  contact_obs$class_dist <-cut(contact_obs$origin_dist, breaks=seq(from=0, to=maxDistance, by=50000), include.lowest = T)
  contact_simul$class_dist <-cut(contact_simul$origin_dist, breaks=seq(from=0, to=maxDistance, by=50000), include.lowest = T)
  
  # Select only contacts with no duplicated enhancers
  obs_stats <- enhancer.statistics[[sp]][[enh]][["real"]]
  simul_stats <- enhancer.statistics[[sp]][[enh]][["simulated"]]
  obs_stats$enh <-  do.call(paste,c(obs_stats[c("chr","start","end")],sep=":"))
  simul_stats$enh <-  do.call(paste,c(simul_stats[c("chr","start","end")],sep=":"))
  obs_stats <- obs_stats[which(obs_stats$BLAT_match < 2),]
  simul_stats <- simul_stats[which(simul_stats$BLAT_match < 2),]
  
  all_obs <- all_obs[which(all_obs$enhancer %in% obs_stats$enh),]
  all_simul <- all_simul[which(all_simul$enhancer %in% simul_stats$enh),]
  contact_obs <- contact_obs[which(contact_obs$origin_enh %in% obs_stats$enh),]
  contact_simul <- contact_simul[which(contact_simul$origin_enh %in% simul_stats$enh),]
  
  
  # Calculate proportion
  conserv <- data.frame(result = sapply(levels(all_obs$class_dist), function(x)
    (nrow(contact_obs[which(contact_obs$class_dist == x ),])/nrow(all_obs[which(all_obs$class_dist == x ),]))*100))
  
  conserv$obs_conflow <- sapply(levels(all_obs$class_dist), function(x)  
    (prop.test(x = nrow(contact_obs[which(contact_obs$class_dist == x ),]),
               n=nrow(all_obs[which(all_obs$class_dist == x ),]), p=0.5)$conf.int[1])*100)
  
  conserv$obs_confup <- sapply(levels(all_obs$class_dist), function(x)  
    (prop.test(x = nrow(contact_obs[which(contact_obs$class_dist == x ),]),
               n=nrow(all_obs[which(all_obs$class_dist == x ),]), p=0.5)$conf.int[2])*100)
 
  conserv$simul <- sapply(levels(all_simul$class_dist), function(x)
    (nrow(contact_simul[which(contact_simul$class_dist == x ),])/nrow(all_simul[which(all_simul$class_dist == x ),]))*100)
  
  conserv$simul_conflow <- sapply(levels(all_simul$class_dist), function(x)  
    (prop.test(x = nrow(contact_simul[which(contact_simul$class_dist == x ),]),
               n=nrow(all_simul[which(all_simul$class_dist == x ),]), p=0.5)$conf.int[1])*100)
  
  conserv$simul_confup <- sapply(levels(all_simul$class_dist), function(x)  
    (prop.test(x = nrow(contact_simul[which(contact_simul$class_dist == x ),]),
               n=nrow(all_simul[which(all_simul$class_dist == x ),]), p=0.5)$conf.int[2])*100)
  
  assign(paste("conserv_dist", enh, sep = "_"), conserv)

}


################################################# Conserv ~ nb cell #################################################

for (enh in enhancer.datasets[[sp]]){
  message(enh)
  all_obs <- gene.enhancer.contacts[[sp]][[enh]][["real"]] 
  all_simul <- gene.enhancer.contacts[[sp]][[enh]][["simulated"]] 
  
  contact_obs <- contacts.conservation[[sp]][[enh]][["real"]] 
  contact_simul <- contacts.conservation[[sp]][[enh]][["simulated"]] 
  
  # Select specific genes
  if(selected_genes != "all"){
    all_obs <- all_obs[which(all_obs$gene %get% dev_gene$Gene.ID),]
    all_simul <- all_simul[which(all_simul$gene %get% dev_gene$Gene.ID),]
    contact_obs <- contact_obs[which(contact_obs$origin_gene %get% dev_gene$Gene.ID),]
    contact_simul <- contact_simul[which(contact_simul$origin_gene %get% dev_gene$Gene.ID),]
  }
  
  
  # Select only contacts with no duplicated enhancers
  obs_stats <- enhancer.statistics[[sp]][[enh]][["real"]]
  simul_stats <- enhancer.statistics[[sp]][[enh]][["simulated"]]
  obs_stats$enh <-  do.call(paste,c(obs_stats[c("chr","start","end")],sep=":"))
  simul_stats$enh <-  do.call(paste,c(simul_stats[c("chr","start","end")],sep=":"))
  obs_stats <- obs_stats[which(obs_stats$BLAT_match < 2),]
  simul_stats <- simul_stats[which(simul_stats$BLAT_match < 2),]
  
  all_obs <- all_obs[which(all_obs$enhancer %in% obs_stats$enh),]
  all_simul <- all_simul[which(all_simul$enhancer %in% simul_stats$enh),]
  contact_obs <- contact_obs[which(contact_obs$origin_enh %in% obs_stats$enh),]
  contact_simul <- contact_simul[which(contact_simul$origin_enh %in% simul_stats$enh),]
  
  # Sample classes
  if(sp=="human"){
    sample_class = 0:7
    class_lab=c(as.character(1:7), ">7")
  }else{sample_class = 0:5
  class_lab=c(as.character(1:5), ">5")}
  
  all_obs$class_dist <- cut(all_obs$nb_cell, breaks=c(sample_class, max(all_obs$nb_cell)), include.lowest = T)
  levels(all_obs$class_dist) = class_lab
  all_simul$class_dist <- cut(all_simul$nb_cell, breaks=c(sample_class, max(all_simul$nb_cell)), include.lowest = T)
  levels(all_simul$class_dist) = class_lab
  
  contact_obs$class_dist <- cut(contact_obs$nb_cell, breaks=c(sample_class, max(contact_obs$nb_cell)), include.lowest = T)
  levels(contact_obs$class_dist) = class_lab
  contact_simul$class_dist <- cut(contact_simul$nb_cell, breaks=c(sample_class, max(contact_simul$nb_cell)), include.lowest = T)
  levels(contact_simul$class_dist) = class_lab
  
  # Overlap with target enhancer
  if (enh %in% c("FANTOM5", "ENCODE")){
    overlap_enh <- read.table(paste(path_overlap, enh, "/enhancer_overlap_target_enhancer.txt", sep=""), header=T, sep="\t")
    enh_tot = nrow(overlap_enh)
    overlap_enh <- overlap_enh[which(overlap_enh[,5] != "NA"),]
    enh_overlap = nrow(overlap_enh)
    message("Proportion of lifted_enh overlap target enh : ", enh_overlap, " on ", enh_tot, " = ", enh_overlap/enh_tot  )
    
    all_obs <- all_obs[which(all_obs$enhancer %in% overlap_enh[,1]),]
    all_simul <- all_simul[which(all_simul$enhancer %in% overlap_enh[,1]),]
    contact_obs <- contact_obs[which(contact_obs$origin_enh %in% overlap_enh[,1]),]
    contact_simul <- contact_simul[which(contact_simul$origin_enh %in% overlap_enh[,1]),]
  }
  
  # Calculate proportion
  conserv_obs <- data.frame(result = sapply(levels(all_obs$class_dist), function(x)
    (nrow(contact_obs[which(contact_obs$class_dist == x),])/(nrow(all_obs[which(all_obs$class_dist == x ),])))*100))
  
  conserv_obs$conflow <- sapply(levels(all_obs$class_dist), function(x)  
    (prop.test(x = nrow(contact_obs[which(contact_obs$class_dist == x ),]),
               n=nrow(all_obs[which(all_obs$class_dist == x ),]), p=0.5)$conf.int[1])*100)
  
  conserv_obs$confup <- sapply(levels(all_obs$class_dist), function(x)  
    (prop.test(x = nrow(contact_obs[which(contact_obs$class_dist == x ),]),
               n=nrow(all_obs[which(all_obs$class_dist == x ),]), p=0.5)$conf.int[2])*100)
  
  conserv_obs$data <- "obs"
  conserv_obs$class <-  class_lab

  conserv_obs$count <- sapply(levels(all_obs$class_dist), function(x) paste("n=", nrow(all_obs[which(all_obs$class_dist == x ),]), sep=""))
  
  conserv_simul <- data.frame(result = sapply(levels(all_simul$class_dist), function(x)
    (nrow(contact_simul[which(contact_simul$class_dist == x ),])/(nrow(all_simul[which(all_simul$class_dist == x ),])))*100))
  
  conserv_simul$conflow <- sapply(levels(all_simul$class_dist), function(x) 
    (prop.test(x = nrow(contact_simul[which(contact_simul$class_dist == x ),]),
               n=nrow(all_simul[which(all_simul$class_dist == x ),]), p=0.5)$conf.int[1])*100)

  conserv_simul$confup <- sapply(levels(all_simul$class_dist), function(x)
    (prop.test(x = nrow(contact_simul[which(contact_simul$class_dist == x ),]), 
               n=nrow(all_simul[which(all_simul$class_dist == x ),]), p=0.5)$conf.int[2])*100)
  
  conserv_simul$data <- "simul"
  conserv_simul$class <- class_lab

  conserv_simul$count <- sapply(levels(all_simul$class_dist), function(x) paste("n=", nrow(all_simul[which(all_simul$class_dist == x ),]), sep=""))
  conserv <- rbind(conserv_obs, conserv_simul)
  conserv$class <- ordered(conserv$class, levels = class_lab)
  
  assign(paste("conserv_sample", enh, sep = "_"), conserv)
  
}

################################################# Conserv in similar samples #################################################
for (enh in enhancer.datasets[[sp]]){
  message(enh)
  data <- c()
  conf_low <- c()
  conf_up <- c()
  p_test <- c()
  n_total <- c()
  
  all_obs <- gene.enhancer.contacts[[sp]][[enh]][["real"]] 
  all_simul <- gene.enhancer.contacts[[sp]][[enh]][["simulated"]] 
  
  contact_obs <- contacts.conservation[[sp]][[enh]][["real"]] 
  contact_simul <- contacts.conservation[[sp]][[enh]][["simulated"]] 
  
  # Select specific genes
  if(selected_genes != "all"){
    all_obs <- all_obs[which(all_obs$gene %get% dev_gene$Gene.ID),]
    all_simul <- all_simul[which(all_simul$gene %get% dev_gene$Gene.ID),]
    contact_obs <- contact_obs[which(contact_obs$origin_gene %get% dev_gene$Gene.ID),]
    contact_simul <- contact_simul[which(contact_simul$origin_gene %get% dev_gene$Gene.ID),]
  }
  
  # Select only contacts with no duplicated enhancers
  obs_stats <- enhancer.statistics[[sp]][[enh]][["real"]]
  simul_stats <- enhancer.statistics[[sp]][[enh]][["simulated"]]
  obs_stats$enh <-  do.call(paste,c(obs_stats[c("chr","start","end")],sep=":"))
  simul_stats$enh <-  do.call(paste,c(simul_stats[c("chr","start","end")],sep=":"))
  obs_stats <- obs_stats[which(obs_stats$BLAT_match < 2),]
  simul_stats <- simul_stats[which(simul_stats$BLAT_match < 2),]
  
  all_obs <- all_obs[which(all_obs$enhancer %in% obs_stats$enh),]
  all_simul <- all_simul[which(all_simul$enhancer %in% simul_stats$enh),]
  contact_obs <- contact_obs[which(contact_obs$origin_enh %in% obs_stats$enh),]
  contact_simul <- contact_simul[which(contact_simul$origin_enh %in% simul_stats$enh),]
  
  # Overlap with target enhancer
  if (enh %in% c("FANTOM5", "ENCODE")){
    overlap_enh <- read.table(paste(path_overlap, enh, "/enhancer_overlap_target_enhancer.txt", sep=""), header=T, sep="\t")
    enh_tot = nrow(overlap_enh)
    overlap_enh <- overlap_enh[which(overlap_enh[,5] != "NA"),]
    enh_overlap = nrow(overlap_enh)
    message("Proportion of lifted_enh overlap target enh : ", enh_overlap, " on ", enh_tot, " = ", enh_overlap/enh_tot  )
    
    all_obs <- all_obs[which(all_obs$enhancer %in% overlap_enh[,1]),]
    all_simul <- all_simul[which(all_simul$enhancer %in% overlap_enh[,1]),]
    contact_obs <- contact_obs[which(contact_obs$origin_enh %in% overlap_enh[,1]),]
    contact_simul <- contact_simul[which(contact_simul$origin_enh %in% overlap_enh[,1]),]
  }
  
  
  ### Select relative cell type
  if (sp == "human"){
    # Pre-adipocytes
    all_obs_adip <- all_obs[which(!is.na(all_obs$pre_adipo)),]
    all_simul_adip <- all_simul[which(!is.na(all_simul$pre_adipo)),]

    # ESC
    all_obs_ESC <- all_obs[which(!is.na(all_obs$hESC)),]
    all_simul_ESC <- all_simul[which(!is.na(all_simul$hESC)),]

    # Bcell
    all_obs_Bcell <- all_obs[which(!is.na(all_obs$TB) | !is.na(all_obs$NB)),]
    all_simul_Bcell <- all_simul[which(!is.na(all_simul$Bcell) | !is.na(all_simul$TB)),]

  }else{
    # Pre-adipocytes
    all_obs_adip <- all_obs[which(!is.na(all_obs$preadip_D0) | !is.na(all_obs$preadip_D2) | !is.na(all_obs$preadip_4H) ),]
    all_simul_adip <-  all_simul[which(!is.na(all_simul$preadip_D0) | !is.na(all_simul$preadip_D2) | !is.na(all_simul$preadip_4H) ),]

    # ESC
    all_obs_ESC <- all_obs[which(!is.na(all_obs$ESC) | !is.na(all_obs$ESC_18)  | !is.na(all_obs$ESC_wild)),]
    all_simul_ESC <- all_simul[which(!is.na(all_simul$ESC) | !is.na(all_simul$ESC_18) | !is.na(all_simul$ESC_wild)),]

    # Bcell
    all_obs_Bcell <- all_obs[which(!is.na(all_obs$preB_aged) | !is.na(all_obs$preB_young)),]
    all_simul_Bcell <- all_simul[which(!is.na(all_simul$preB_aged) | !is.na(all_simul$preB_young)),]
    }
  
  # Pre-adipocytes
  contact_obs_adip <- contact_obs[which(contact_obs$pre_adipo > 0 & (contact_obs$preadip_D0 > 0 | contact_obs$preadip_D2 > 0 | contact_obs$preadip_4H > 0 )),]
  contact_simul_adip <- contact_simul[which(contact_simul$pre_adipo > 0 & (contact_simul$preadip_D0 > 0 | contact_simul$preadip_D2 > 0 | contact_simul$preadip_4H > 0 )),]
  
  # ESC
  contact_obs_ESC <- contact_obs[which(contact_obs$hESC > 0 & (contact_obs$ESC > 0 | contact_obs$ESC_18 > 0  | contact_obs$ESC_wild > 0 )),]
  contact_simul_ESC <- contact_simul[which(contact_simul$hESC > 0 & (contact_simul$ESC > 0 | contact_simul$ESC_18 > 0 | contact_simul$ESC_wild > 0)),]
  
  # Bcell
  contact_obs_Bcell <- contact_obs[which((contact_obs$TB > 0 | contact_obs$NB > 0) & (contact_obs$preB_aged > 0 | contact_obs$preB_young > 0)),]
  contact_simul_Bcell <- contact_simul[which((contact_simul$TB > 0 | contact_simul$NB > 0) & (contact_simul$preB_aged > 0 | contact_simul$preB_young > 0)),]
  
  for (cell in c("adip", "ESC", "Bcell")){
    # Calculate proportion
    all_obs = get(paste("all_obs", cell, sep="_"))
    all_simul = get(paste("all_simul", cell, sep="_"))
    contact_obs = get(paste("contact_obs", cell, sep="_"))
    contact_simul = get(paste("contact_simul", cell, sep="_"))
    
    mat <- matrix(c(nrow(contact_obs), nrow(contact_simul), nrow(all_obs)-nrow(contact_obs), nrow(all_simul)-nrow(contact_simul)),2)
    
    conf_low <- append(conf_low, c((prop.test(x = nrow(contact_obs), n=nrow(all_obs), p=0.5)$conf.int[1])*100,
                                   (prop.test(x = nrow(contact_simul), n=nrow(all_simul), p=0.5)$conf.int[1])*100, NA))
    
    conf_up <- append(conf_up, c((prop.test(x = nrow(contact_obs), n=nrow(all_obs), p=0.5)$conf.int[2])*100,
                                 (prop.test(x = nrow(contact_simul), n=nrow(all_simul), p=0.5)$conf.int[2])*100, NA))
    
    data <- append(data, c( (nrow(contact_obs)/nrow(all_obs))*100, (nrow(contact_simul)/nrow(all_simul))*100, NA))
    
    n_total <- append(n_total, c(paste0("N = ", nrow(all_obs)), paste0("N = ", nrow(all_simul)), NA))
    
  }

  conserv <- data.frame(data=data, conf_low=conf_low, conf_up=conf_up, n_total=n_total)
  assign(paste("conserv_similar_sample", enh, sep = "_"), conserv)
  
}

if (sp == "mouse"){
  save(conserv_global,
       conserv_dist_FANTOM5, conserv_dist_ENCODE, 
       conserv_sample_FANTOM5, conserv_sample_ENCODE,
       conserv_similar_sample_FANTOM5, conserv_similar_sample_ENCODE,
       file = paste(pathFigures, "RData/Fig5_", sp, "_", selected_genes, "_genes.Rdata", sep=""))
  
}else{
  save(conserv_global,
       conserv_dist_FANTOM5, conserv_dist_ENCODE, conserv_dist_RoadmapEpigenomics, conserv_dist_FOCS_GRO_seq,
       conserv_sample_FANTOM5, conserv_sample_ENCODE, conserv_sample_RoadmapEpigenomics, conserv_sample_FOCS_GRO_seq,
       conserv_similar_sample_FANTOM5, conserv_similar_sample_ENCODE, conserv_similar_sample_RoadmapEpigenomics, conserv_similar_sample_FOCS_GRO_seq,
       file = paste(pathFigures, "RData/Fig5_", sp, "_", selected_genes, "_genes.Rdata", sep=""))
}


