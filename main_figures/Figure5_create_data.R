################################################################################################################################################
library(data.table)
options(stringsAsFactors = FALSE)

source("parameters.R") ## pathFinalDataFinalDatas are defined based on the user name

# Choose species within : human ; mouse
ref_sp = "human" 

# Choose genes within : all ; dvpt ; other
selected_genes = "all"
if (selected_genes == "dvpt"){`%get%` <- `%in%`}else{`%get%` <- Negate(`%in%`)}

minDistance=25e3
maxDistance=2e6

enhancers = c("FANTOM5", "ENCODE")

if(ref_sp == "human"){
  enhancers <- c(enhancers, "RoadmapEpigenomics", "FOCS_GRO_seq")
  target_sp = "mouse"
}

path_evol <- paste(pathFinalData, "SupplementaryDataset7/", ref_sp, "/", sep="")
path_contact <- paste(pathFinalData, "SupplementaryDataset4/", ref_sp, "/", sep="")


##############################################################################################################################
################################################# Conserved contact global   #################################################
data_global <- c()
conf_low_global <- c()
conf_up_global <- c()
n_total_global <- c()

for (enh in enhancers){
  message(enh)
  all_obs <- read.table(paste(path_contact, "/", enh, "/gene_enhancer_contacts_original_interactions.txt", sep=""), header=T, sep="\t")
  all_simul <- read.table(paste(path_contact, "/", enh, "/gene_enhancer_contacts_simulated_interactions.txt", sep=""), header=T, sep="\t")
  contact_obs <- read.table(paste(path_evol,"/contact_conservation/", enh, "/", ref_sp, "_original2", target_sp, "_original.txt", sep=""), header=T, sep="\t")
  contact_simul <- read.table(paste(path_evol,"/contact_conservation/", enh, "/", ref_sp, "_simulated2", target_sp, "_simulated.txt", sep=""), header=T, sep="\t")
  
  # Select specific genes
  if(selected_genes != "all"){
    dev_gene <- read.table(paste(pathFinalData, "SupplementaryDataset3/genes/", ref_sp, "_dvpt_process_genes.txt", sep=""), header=T, sep="\t")
    all_obs <- all_obs[which(all_obs$gene %get% dev_gene$Gene.ID),]
    all_simul <- all_simul[which(all_simul$gene %get% dev_gene$Gene.ID),]
    contact_obs <- contact_obs[which(contact_obs$origin_gene %get% dev_gene$Gene.ID),]
    contact_simul <- contact_simul[which(contact_simul$origin_gene %get% dev_gene$Gene.ID),]
  }
  
  
  # Select only contacts with no duplicated enhancers
  obs_stats <- read.table(paste(path_contact, enh, "/statistics_contacted_enhancers_original.txt", sep=""), header=T)
  simul_stats <- read.table(paste(path_contact, enh, "/statistics_contacted_enhancers_simulated.txt", sep=""), header=T)
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
    overlap_enh <- read.table(paste(path_contact, enh, "/enhancer_overlap_target_enhancer.txt", sep=""), header=T, sep="\t")
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
  mat <- matrix(c(nrow(contact_obs), nrow(contact_simul), nrow(all_obs)-nrow(contact_obs), nrow(all_simul)-nrow(contact_simul)),2)
  
  conf_low_global <- append(conf_low_global, c((prop.test(x = nrow(contact_obs), n=nrow(all_obs), p=0.5)$conf.int[1])*100,
                                 (prop.test(x = nrow(contact_simul), n=nrow(all_simul), p=0.5)$conf.int[1])*100, NA))
  
  conf_up_global <- append(conf_up_global, c((prop.test(x = nrow(contact_obs), n=nrow(all_obs), p=0.5)$conf.int[2])*100,
                               (prop.test(x = nrow(contact_simul), n=nrow(all_simul), p=0.5)$conf.int[2])*100, NA))
  
  data_global <- append(data_global, c( (nrow(contact_obs)/nrow(all_obs))*100, (nrow(contact_simul)/nrow(all_simul))*100, NA))
  
  n_total_global <- append(n_total_global, c(paste0("N = ", nrow(all_obs)), paste0("N = ", nrow(all_simul)), NA))
  
  conserv <- data.frame(data=data_global, conf_low=conf_low_global, conf_up=conf_up_global, n_total=n_total_global)
  assign("conserv_global", conserv)
  
}

################################################# Conserv ~ genomic distance   #################################################

for (enh in enhancers){
  all_obs <- read.table(paste(path_contact, "/", enh, "/gene_enhancer_contacts_original_interactions.txt", sep=""), header=T, sep="\t")
  all_simul <- read.table(paste(path_contact, "/", enh, "/gene_enhancer_contacts_simulated_interactions.txt", sep=""), header=T, sep="\t")
  contact_obs <- read.table(paste(path_evol,"/contact_conservation/", enh, "/", ref_sp, "_original2", target_sp, "_original.txt", sep=""), header=T, sep="\t")
  contact_simul <- read.table(paste(path_evol,"/contact_conservation/", enh, "/", ref_sp, "_simulated2", target_sp, "_simulated.txt", sep=""), header=T, sep="\t")
  
  # Select specific genes
  if(selected_genes != "all"){
    dev_gene <- read.table(paste(pathFinalData, "SupplementaryDataset3/genes/", ref_sp, "_dvpt_process_genes.txt", sep=""), header=T, sep="\t")
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
  obs_stats <- read.table(paste(path_contact, enh, "/statistics_contacted_enhancers_original.txt", sep=""), header=T)
  simul_stats <- read.table(paste(path_contact, enh, "/statistics_contacted_enhancers_simulated.txt", sep=""), header=T)
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


################################################# Conserv ~ nb samples #################################################

for (enh in enhancers){
  message(enh)
  all_obs <- read.table(paste(path_contact, "/", enh, "/gene_enhancer_contacts_original_interactions.txt", sep=""), header=T, sep="\t")
  all_simul <- read.table(paste(path_contact, "/", enh, "/gene_enhancer_contacts_simulated_interactions.txt", sep=""), header=T, sep="\t")
  contact_obs <- read.table(paste(path_evol,"/contact_conservation/", enh, "/", ref_sp, "_original2", target_sp, "_original.txt", sep=""), header=T, sep="\t")
  contact_simul <- read.table(paste(path_evol,"/contact_conservation/", enh, "/", ref_sp, "_simulated2", target_sp, "_simulated.txt", sep=""), header=T, sep="\t")
  
  # Select specific genes
  if(selected_genes != "all"){
    dev_gene <- read.table(paste(pathFinalData, "SupplementaryDataset3/genes/", ref_sp, "_dvpt_process_genes.txt", sep=""), header=T, sep="\t")
    all_obs <- all_obs[which(all_obs$gene %get% dev_gene$Gene.ID),]
    all_simul <- all_simul[which(all_simul$gene %get% dev_gene$Gene.ID),]
    contact_obs <- contact_obs[which(contact_obs$origin_gene %get% dev_gene$Gene.ID),]
    contact_simul <- contact_simul[which(contact_simul$origin_gene %get% dev_gene$Gene.ID),]
  }
  
  # Select only contacts with no duplicated enhancers
  obs_stats <- read.table(paste(path_contact, enh, "/statistics_contacted_enhancers_original.txt", sep=""), header=T)
  simul_stats <- read.table(paste(path_contact, enh, "/statistics_contacted_enhancers_simulated.txt", sep=""), header=T)
  obs_stats$enh <-  do.call(paste,c(obs_stats[c("chr","start","end")],sep=":"))
  simul_stats$enh <-  do.call(paste,c(simul_stats[c("chr","start","end")],sep=":"))
  obs_stats <- obs_stats[which(obs_stats$BLAT_match < 2),]
  simul_stats <- simul_stats[which(simul_stats$BLAT_match < 2),]
  
  all_obs <- all_obs[which(all_obs$enhancer %in% obs_stats$enh),]
  all_simul <- all_simul[which(all_simul$enhancer %in% simul_stats$enh),]
  contact_obs <- contact_obs[which(contact_obs$origin_enh %in% obs_stats$enh),]
  contact_simul <- contact_simul[which(contact_simul$origin_enh %in% simul_stats$enh),]
  
  # Sample classes
  if(ref_sp=="human"){
    sample_class = c(1, 2, 5, 10, 20)
    class_lab=c("1", "2-5", "6-10", "11-20", "21-26")
  }else{sample_class = c(1, 2, 3, 6, 9)
  class_lab=c("1", "2-3", "4-6", "7-9", "10-14")}
  
  all_obs$class_dist <- cut(all_obs$nb_sample, breaks=c(sample_class, max(all_obs$nb_sample)), include.lowest = T)
  all_simul$class_dist <- cut(all_simul$nb_sample, breaks=c(sample_class, max(all_simul$nb_sample)+1), include.lowest = T)
  contact_obs$class_dist <- cut(contact_obs$nb_sample, breaks=c(sample_class, max(contact_obs$nb_sample)), include.lowest = T)
  contact_simul$class_dist <- cut(contact_simul$nb_sample, breaks=c(sample_class, max(contact_simul$nb_sample)), include.lowest = T)
  
  # Overlap with target enhancer
  if (enh %in% c("FANTOM5", "ENCODE")){
    overlap_enh <- read.table(paste(path_contact, enh, "/enhancer_overlap_target_enhancer.txt", sep=""), header=T, sep="\t")
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
    tryCatch((prop.test(x = nrow(contact_simul[which(contact_simul$class_dist == x ),]),
               n=nrow(all_simul[which(all_simul$class_dist == x ),]), p=0.5)$conf.int[1])*100,error=function(e) 0))

  conserv_simul$confup <- sapply(levels(all_simul$class_dist), function(x)  
    tryCatch((prop.test(x = nrow(contact_simul[which(contact_simul$class_dist == x ),]),
                        n=nrow(all_simul[which(all_simul$class_dist == x ),]), p=0.5)$conf.int[2])*100,error=function(e) 0))
  
  conserv_simul$data <- "simul"
  conserv_simul$class <- class_lab

  conserv_simul$count <- sapply(levels(all_simul$class_dist), function(x) paste("n=", nrow(all_simul[which(all_simul$class_dist == x ),]), sep=""))
  conserv <- rbind(conserv_obs, conserv_simul)
  conserv$class <- ordered(conserv$class, levels = class_lab)
  
  assign(paste("conserv_sample", enh, sep = "_"), conserv)
  
}

################################################# Conserv in similar samples #################################################
for (enh in enhancers){
  message(enh)
  data <- c()
  conf_low <- c()
  conf_up <- c()
  p_test <- c()
  n_total <- c()
  
  all_obs <- read.table(paste(path_contact, "/", enh, "/gene_enhancer_contacts_original_interactions.txt", sep=""), header=T, sep="\t")
  all_simul <- read.table(paste(path_contact, "/", enh, "/gene_enhancer_contacts_simulated_interactions.txt", sep=""), header=T, sep="\t")
  contact_obs <- read.table(paste(path_evol,"/contact_conservation/", enh, "/", ref_sp, "_original2", target_sp, "_original.txt", sep=""), header=T, sep="\t")
  contact_simul <- read.table(paste(path_evol,"/contact_conservation/", enh, "/", ref_sp, "_simulated2", target_sp, "_simulated.txt", sep=""), header=T, sep="\t")
  
  # Select specific genes
  if(selected_genes != "all"){
    dev_gene <- read.table(paste(pathFinalData, "SupplementaryDataset3/genes/", ref_sp, "_dvpt_process_genes.txt", sep=""), header=T, sep="\t")
    all_obs <- all_obs[which(all_obs$gene %get% dev_gene$Gene.ID),]
    all_simul <- all_simul[which(all_simul$gene %get% dev_gene$Gene.ID),]
    contact_obs <- contact_obs[which(contact_obs$origin_gene %get% dev_gene$Gene.ID),]
    contact_simul <- contact_simul[which(contact_simul$origin_gene %get% dev_gene$Gene.ID),]
  }
  
  # Select only contacts with no duplicated enhancers
  obs_stats <- read.table(paste(path_contact, enh, "/statistics_contacted_enhancers_original.txt", sep=""), header=T)
  simul_stats <- read.table(paste(path_contact, enh, "/statistics_contacted_enhancers_simulated.txt", sep=""), header=T)
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
    overlap_enh <- read.table(paste(path_contact, enh, "/enhancer_overlap_target_enhancer.txt", sep=""), header=T, sep="\t")
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
  if (ref_sp == "human"){
    # Pre-adipocytes
    all_obs_adip <- all_obs[which(!is.na(all_obs$pre_adipo)),]
    all_simul_adip <- all_simul[which(!is.na(all_simul$pre_adipo)),]

    # ESC
    all_obs_ESC <- all_obs[which(!is.na(all_obs$hESC)),]
    all_simul_ESC <- all_simul[which(!is.na(all_simul$hESC)),]

    # Bcell
    all_obs_Bcell <- all_obs[which(!is.na(all_obs$Bcell) | !is.na(all_obs$TB) | !is.na(all_obs$NB)),]
    all_simul_Bcell <- all_simul[which(!is.na(all_simul$Bcell) | !is.na(all_simul$TB) | !is.na(all_simul$NB)),]

  }else{
    # Pre-adipocytes
    all_obs_adip <- all_obs[which(!is.na(all_obs$preadip_D0) | !is.na(all_obs$preadip_D2) | !is.na(all_obs$preadip_4H) ),]
    all_simul_adip <-  all_simul[which(!is.na(all_simul$preadip_D0) | !is.na(all_simul$preadip_D2) | !is.na(all_simul$preadip_4H) ),]

    # ESC
    all_obs_ESC <- all_obs[which(!is.na(all_obs$ESC) | !is.na(all_obs$ESC_18) | !is.na(all_obs$ESC_NKO) | !is.na(all_obs$ESC_wild)),]
    all_simul_ESC <- all_simul[which(!is.na(all_simul$ESC) | !is.na(all_simul$ESC_18) | !is.na(all_simul$ESC_NKO) | !is.na(all_simul$ESC_wild)),]

    # Bcell
    all_obs_Bcell <- all_obs[which(!is.na(all_obs$preB_aged) | !is.na(all_obs$preB_young)),]
    all_simul_Bcell <- all_simul[which(!is.na(all_simul$preB_aged) | !is.na(all_simul$preB_young)),]
    }
  
  # Pre-adipocytes
  contact_obs_adip <- contact_obs[which(contact_obs$pre_adipo > 0 & (contact_obs$preadip_D0 > 0 | contact_obs$preadip_D2 > 0 | contact_obs$preadip_4H > 0 )),]
  contact_simul_adip <- contact_simul[which(contact_simul$pre_adipo > 0 & (contact_simul$preadip_D0 > 0 | contact_simul$preadip_D2 > 0 | contact_simul$preadip_4H > 0 )),]
  
  # ESC
  contact_obs_ESC <- contact_obs[which(contact_obs$hESC > 0 & (contact_obs$ESC > 0 | contact_obs$ESC_18 > 0 | contact_obs$ESC_NKO > 0 | contact_obs$ESC_wild > 0 )),]
  contact_simul_ESC <- contact_simul[which(contact_simul$hESC > 0 & (contact_simul$ESC > 0 | contact_simul$ESC_18 > 0 | contact_simul$ESC_NKO > 0 | contact_simul$ESC_wild > 0)),]
  
  # Bcell
  contact_obs_Bcell <- contact_obs[which((contact_obs$Bcell > 0 | contact_obs$TB > 0) & (contact_obs$preB_aged > 0 | contact_obs$preB_young > 0)),]
  contact_simul_Bcell <- contact_simul[which((contact_simul$Bcell > 0 | contact_simul$TB > 0) & (contact_simul$preB_aged > 0 | contact_simul$preB_young > 0)),]
  
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

if (ref_sp == "mouse"){
  save(conserv_global,
       conserv_dist_FANTOM5, conserv_dist_ENCODE, 
       conserv_sample_FANTOM5, conserv_sample_ENCODE,
       conserv_similar_sample_FANTOM5, conserv_similar_sample_ENCODE,
       file = paste(pathFigures, "/Fig5_", ref_sp, "_", selected_genes, "_genes_unique.Rdata", sep=""))
  
}else{
  save(conserv_global,
       conserv_dist_FANTOM5, conserv_dist_ENCODE, conserv_dist_RoadmapEpigenomics, conserv_dist_FOCS_GRO_seq,
       conserv_sample_FANTOM5, conserv_sample_ENCODE, conserv_sample_RoadmapEpigenomics, conserv_sample_FOCS_GRO_seq,
       conserv_similar_sample_FANTOM5, conserv_similar_sample_ENCODE, conserv_similar_sample_RoadmapEpigenomics, conserv_similar_sample_FOCS_GRO_seq,
       file = paste(pathFigures, "/Fig5_", ref_sp, "_", selected_genes, "_genes_unique.Rdata", sep=""))
}





# test_obs <- data.frame(result = sapply(levels(as.factor(all_obs$nb_sample)), function(x)
#   (nrow(contact_obs[which(as.factor(contact_obs$nb_sample) == x),])/(nrow(all_obs[which(as.factor(all_obs$nb_sample) == x ),])))*100))
# 
# test_obs$data <- "obs"
# test_obs$class <- levels(as.factor(all_obs$nb_sample))
# 
# test_simul <- data.frame(result = sapply(levels(as.factor(all_simul$nb_sample)), function(x)
#   (nrow(contact_simul[which(as.factor(contact_simul$nb_sample) == x),])/(nrow(all_simul[which(as.factor(all_simul$nb_sample) == x ),])))*100))
# 
# test_simul$data <- "simul"
# test_simul$class <- levels(as.factor(all_simul$nb_sample))
# test_conserv <- rbind(test_obs, test_simul)
# test_conserv$class <- ordered(test_conserv$class, levels = levels(as.factor(all_obs$nb_sample)))
# 
# bar <- barplot(result ~ data+class, beside=T, data=test_conserv, space = c(0.1, 0.6),
#                border=c("darkgreen", "firebrick1"), col="white", ylab="Conserved contact (%)", xlab="Number of sample", 
#                main="")
