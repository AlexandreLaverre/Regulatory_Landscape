################################################################################################################################################
library(naniar)
options(stringsAsFactors = FALSE)

ref_sp = "human" # to change human or mouse
target_sp = "mouse"

source("parameters.R") 
path_evol <- paste(pathFinalData, "SupplementaryDataset7/", ref_sp, "/", sep="")
path_annot <- paste(pathFinalData, "SupplementaryDataset4/", ref_sp, "/", sep="")


minDistance=25e3
maxDistance=2e6

enhancers = c("FANTOM5", "ENCODE")
if(ref_sp == "human"){enhancers <- c(enhancers, "RoadmapEpigenomics", "FOCS_GRO_seq")}
species <- c("macaque", "dog", "cow", "elephant", "rabbit", "rat", target_sp, "opossum", "chicken")

###################################################################################################################################
################################################# Conserv synteny ~ all species   #################################################
conserv_synteny <- list()
for (enh in enhancers){
  data <- c()
  p_test <- c()
  conf_low <- c()
  conf_up <- c()
  message("Running conserv synteny with ", enh)
  for (sp in species){
    synt_obs <- read.table(paste(path_evol,"/synteny_conservation/", enh, "/", ref_sp, "2", sp, "_original_synteny.txt", sep=""), header=T)
    synt_simul <- read.table(paste(path_evol,"/synteny_conservation/", enh, "/", ref_sp, "2", sp, "_simulated_synteny.txt", sep=""), header=T)
    
    # Filters
    synt_obs <- synt_obs[which(synt_obs$BLAT_match < 2 & synt_obs$align_score > 0.4 & synt_obs$origin_dist < maxDistance),]
    synt_simul <- synt_simul[which(synt_simul$BLAT_match < 2 & synt_simul$align_score > 0.4 & synt_simul$origin_dist < maxDistance),]
    synt_obs <- synt_obs %>% replace_with_na(replace = list(target_dist = "trans"))
    synt_simul <- synt_simul %>% replace_with_na(replace = list(target_dist = "trans"))
    
    # Calculate proportion
    nb_synt_obs = nrow(synt_obs[which(as.numeric(synt_obs$target_dist) < maxDistance),])
    nb_synt_simul = nrow(synt_simul[which(as.numeric(synt_simul$target_dist) < maxDistance),])
    
    mat <- matrix(c(nb_synt_obs, nb_synt_simul, nrow(synt_obs)-nb_synt_obs, nrow(synt_simul)-nb_synt_simul),2)
    p_test <- append(p_test, prop.test(mat)$p.value)
    
    conf_low <- append(conf_low, c((prop.test(x = nb_synt_obs, n=nrow(synt_obs), p=0.5)$conf.int[1]),
                                   (prop.test(x = nb_synt_simul, n=nrow(synt_simul), p=0.5)$conf.int[1]), NA))
    
    conf_up <- append(conf_up, c((prop.test(x = nb_synt_obs, n=nrow(synt_obs), p=0.5)$conf.int[2]),
                                 (prop.test(x = nb_synt_simul, n=nrow(synt_simul), p=0.5)$conf.int[2]), NA))
    
    data <- append(data, c(nb_synt_obs/nrow(synt_obs), nb_synt_simul/nrow(synt_simul), NA))
    
    message(sp, " done !")
  }
  
  conserv <- data.frame(data=data, conf_low=conf_low, conf_up=conf_up, p_test=p_test)
  conserv_synteny[[enh]] <- conserv

}  

########################################################################################################################################
################################################# Conserv synteny ~ genomic distance   #################################################
for (enh in enhancers){
  conserv_dist_list <- list()
  message("#### Synteny according to dist for ", enh, " enhancers ####")
  for (sp in species){
    synt_obs <- read.table(paste(path_evol,"/synteny_conservation/", enh, "/", ref_sp, "2", sp, "_original_synteny.txt", sep=""), header=T)
    synt_simul <- read.table(paste(path_evol,"/synteny_conservation/", enh, "/", ref_sp, "2", sp, "_simulated_synteny.txt", sep=""), header=T)
    
    # Filters
    synt_obs <- synt_obs[which(synt_obs$BLAT_match < 2  & synt_obs$align_score > 0.4 & synt_obs$origin_dist < maxDistance),]
    synt_simul <- synt_simul[which(synt_simul$BLAT_match < 2  & synt_simul$align_score > 0.4 & synt_simul$origin_dist < maxDistance),]
    synt_obs <- synt_obs %>% replace_with_na(replace = list(target_dist = "trans"))
    synt_simul <- synt_simul %>% replace_with_na(replace = list(target_dist = "trans"))
    
    # Distance classes
    synt_obs$class_dist <-cut(synt_obs$origin_dist, breaks=seq(from=0, to=maxDistance, by=50000), include.lowest = T)
    synt_simul$class_dist <-cut(synt_simul$origin_dist, breaks=seq(from=0, to=maxDistance, by=50000), include.lowest = T)
    
    ## Calculate proportion
    # Original
    conserv <- data.frame(obs = sapply(levels(synt_obs$class_dist), function(x)
      nrow(synt_obs[which(as.numeric(synt_obs$target_dist) < maxDistance & synt_obs$class_dist == x ),])
      /nrow(synt_obs[which(synt_obs$class_dist == x ),])))
    
    conserv$conf_low_obs <- sapply(levels(synt_obs$class_dist), function(x)
      tryCatch(prop.test(x = nrow(synt_obs[which(as.numeric(synt_obs$target_dist) < maxDistance & synt_obs$class_dist == x ),]),
                n=nrow(synt_obs[which(synt_obs$class_dist == x ),]), p=0.5)$conf.int[1], error=function(e) NA))
    
    conserv$conf_up_obs <-sapply(levels(synt_obs$class_dist), function(x)
      tryCatch(prop.test(x = nrow(synt_obs[which(as.numeric(synt_obs$target_dist) < maxDistance & synt_obs$class_dist == x ),]),
                n=nrow(synt_obs[which(synt_obs$class_dist == x ),]), p=0.5)$conf.int[2], error=function(e) NA))
    
    # Simulated
    conserv$simul <- sapply(levels(synt_simul$class_dist), function(x)
      nrow(synt_simul[which(as.numeric(synt_simul$target_dist) < maxDistance & synt_simul$class_dist == x ),])
      /nrow(synt_simul[which(synt_simul$class_dist == x ),]))
    

    conserv$conf_low_simul <- sapply(levels(synt_simul$class_dist), function(x)
      tryCatch(prop.test(x = nrow(synt_simul[which(as.numeric(synt_simul$target_dist) < maxDistance & synt_simul$class_dist == x ),]),
                n=nrow(synt_simul[which(synt_simul$class_dist == x ),]), p=0.5)$conf.int[1], error=function(e) NA))
    
    conserv$conf_up_simul <- sapply(levels(synt_simul$class_dist), function(x)
      tryCatch(prop.test(x = nrow(synt_simul[which(as.numeric(synt_simul$target_dist) < maxDistance & synt_simul$class_dist == x ),]),
                n=nrow(synt_simul[which(synt_simul$class_dist == x ),]), p=0.5)$conf.int[2], error=function(e) NA))
      
    conserv_dist_list[[sp]] <- conserv
    message(sp, " done !")
  }
  
  assign(paste("conserv_synteny_dist", enh, sep = "_"), conserv_dist_list)

}  

################################################################################################################################################
# Output

if (ref_sp == "mouse"){
  save(conserv_synteny, conserv_synteny_dist_FANTOM5, conserv_synteny_dist_ENCODE, file = paste(pathFigures, "/Fig4_", ref_sp, ".Rdata", sep=""))
}else{
  save(conserv_synteny, conserv_synteny_dist_FANTOM5, conserv_synteny_dist_ENCODE, conserv_synteny_dist_RoadmapEpigenomics, conserv_synteny_dist_FOCS_GRO_seq,
       file = paste(pathFigures, "Rdata/Fig4_", ref_sp, ".Rdata", sep=""))
}