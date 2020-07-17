options(stringsAsFactors = FALSE)

ref_sp = "mouse"
target_sp = "human"

enhancers <- c("CAGE", "ENCODE")
if(ref_sp == "human"){enhancers <- c(enhancers, "RoadMap", "GRO_seq")}

species <- c("macaque", "dog", "cow", "elephant", "rabbit", "rat", target_sp, "opossum", "chicken")

path <- "/home/laverre/Data/Regulatory_landscape/result"
path_evol <- paste(path, "Supplementary_dataset6_regulatory_landscape_evolution", ref_sp, sep="/")
path_annot <- paste(path, "Supplementary_dataset3_annotations", ref_sp, sep="/")

max_dist = 2000000

################################################# Conserv synteny ~ all species   #################################################
library(naniar)

conserv_synteny <- list()
for (enh in enhancers){
  data <- c()
  p_test <- c()
  conf_low <- c()
  conf_up <- c()
  message("Running conserv synteny with ", enh)
  for (sp in species){
    synt_obs <- read.table(paste(path_evol,"/synteny_conservation/", enh, "/", ref_sp, "2", sp, "_", enh, "_original_synteny.txt", sep=""), header=T)
    synt_simul <- read.table(paste(path_evol,"/synteny_conservation/", enh, "/", ref_sp, "2", sp, "_", enh, "_simulated_synteny.txt", sep=""), header=T)
    
    # Filters
    synt_obs <- synt_obs[which(synt_obs$BLAT_match < 2 & synt_obs$align_score > 0.4 & synt_obs$origin_dist < max_dist),]
    synt_simul <- synt_simul[which(synt_simul$BLAT_match < 2 & synt_simul$align_score > 0.4 & synt_simul$origin_dist < max_dist),]
    synt_obs <- synt_obs %>% replace_with_na(replace = list(target_dist = "trans"))
    synt_simul <- synt_simul %>% replace_with_na(replace = list(target_dist = "trans"))
    
    # Calculate proportion
    nb_synt_obs = nrow(synt_obs[which(!is.na(synt_obs$target_dist) & as.numeric(synt_obs$target_dist) < max_dist),])
    nb_synt_simul = nrow(synt_simul[which(!is.na(synt_simul$target_dist) & as.numeric(synt_simul$target_dist) < max_dist),])
    
    nb_synt_obs = nrow(synt_obs[which(as.numeric(synt_obs$target_dist) < max_dist),])
    nb_synt_simul = nrow(synt_simul[which(as.numeric(synt_simul$target_dist) < max_dist),])
    
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


################################################# Conserv synteny ~ genomic distance   #################################################
for (enh in enhancers){
  conserv_dist_list <- list()
  message("#### Synteny according to dist for ", enh, " enhancers ####")
  for (sp in species){
    
    synt_obs <- read.table(paste(path_evol,"/synteny_conservation/", enh, "/", ref_sp, "2", sp, "_", enh, "_original_synteny.txt", sep=""), header=T)
    synt_simul <- read.table(paste(path_evol,"/synteny_conservation/", enh, "/", ref_sp, "2", sp, "_", enh, "_simulated_synteny.txt", sep=""), header=T)
    
    # Filters
    synt_obs <- synt_obs[which(synt_obs$BLAT_match < 2  & synt_obs$align_score > 0.4 & synt_obs$origin_dist < max_dist),]
    synt_simul <- synt_simul[which(synt_simul$BLAT_match < 2  & synt_simul$align_score > 0.4 & synt_simul$origin_dist < max_dist),]
    synt_obs <- synt_obs %>% replace_with_na(replace = list(target_dist = "trans"))
    synt_simul <- synt_simul %>% replace_with_na(replace = list(target_dist = "trans"))
    
    # Distance classes
    synt_obs$class_dist <-cut(synt_obs$origin_dist, breaks=seq(from=0, to=max_dist, by=50000), include.lowest = T)
    synt_simul$class_dist <-cut(synt_simul$origin_dist, breaks=seq(from=0, to=max_dist, by=50000), include.lowest = T)
    class_leg <- c("0", "500Kb", "1Mb", "1.5Mb", "2Mb")
    
    # Calculate proportion
    conserv <- data.frame(obs = sapply(levels(synt_obs$class_dist), function(x)
      nrow(synt_obs[which(!is.na(synt_obs$target_dist) & as.numeric(synt_obs$target_dist) < max_dist & synt_obs$class_dist == x ),])
      /nrow(synt_obs[which(synt_obs$class_dist == x ),])))
    
    conserv$conf_low_obs <- sapply(levels(synt_obs$class_dist), function(x)
      tryCatch(prop.test(x = nrow(synt_obs[which(!is.na(synt_obs$target_dist) & as.numeric(synt_obs$target_dist) < max_dist & synt_obs$class_dist == x ),]),
                n=nrow(synt_obs[which(synt_obs$class_dist == x ),]), p=0.5)$conf.int[1], error=function(e) NA))
    
    conserv$conf_up_obs <-sapply(levels(synt_obs$class_dist), function(x)
      tryCatch(prop.test(x = nrow(synt_obs[which(!is.na(synt_obs$target_dist) & as.numeric(synt_obs$target_dist) < max_dist & synt_obs$class_dist == x ),]),
                n=nrow(synt_obs[which(synt_obs$class_dist == x ),]), p=0.5)$conf.int[2], error=function(e) NA))
    
    conserv$simul <- sapply(levels(synt_simul$class_dist), function(x)
      nrow(synt_simul[which(!is.na(synt_simul$target_dist) & as.numeric(synt_simul$target_dist) < max_dist & synt_simul$class_dist == x ),])
      /nrow(synt_simul[which(synt_simul$class_dist == x ),]))
    

    conserv$conf_low_simul <- sapply(levels(synt_simul$class_dist), function(x)
      tryCatch(prop.test(x = nrow(synt_simul[which(!is.na(synt_simul$target_dist) & as.numeric(synt_simul$target_dist) < max_dist & synt_simul$class_dist == x ),]),
                n=nrow(synt_simul[which(synt_simul$class_dist == x ),]), p=0.5)$conf.int[1], error=function(e) NA))
    
    conserv$conf_up_simul <- sapply(levels(synt_simul$class_dist), function(x)
      tryCatch(prop.test(x = nrow(synt_simul[which(!is.na(synt_simul$target_dist) & as.numeric(synt_simul$target_dist) < max_dist & synt_simul$class_dist == x ),]),
                n=nrow(synt_simul[which(synt_simul$class_dist == x ),]), p=0.5)$conf.int[2], error=function(e) NA))
      
    conserv_dist_list[[sp]] <- conserv
    message(sp, " done !")
  }
  
  assign(paste("conserv_synteny_dist", enh, sep = "_"), conserv_dist_list)

}  

if (ref_sp == "mouse"){
  save(conserv_synteny, conserv_synteny_dist_CAGE, conserv_synteny_dist_ENCODE, file = paste(path, "/Figures/Fig4_", ref_sp, ".Rdata", sep=""))
}else{
  save(conserv_synteny, conserv_synteny_dist_CAGE, conserv_synteny_dist_ENCODE, conserv_synteny_dist_RoadMap, conserv_synteny_dist_GRO_seq,
       file = paste(path, "/Figures/Fig4_", ref_sp, ".Rdata", sep=""))
}


########################## Intersection human - opposum & human - mouse ########################## 
enh='ENCODE'

for (sp in c("opossum", "mouse")){
  synt_obs <- read.table(paste(path_evol,"/synteny_conservation/", enh, "/human2", sp, "_", enh, "_original_synteny.txt", sep=""), header=T)
  synt_simul <- read.table(paste(path_evol,"/synteny_conservation/", enh, "/human2", sp, "_", enh, "_simulated_synteny.txt", sep=""), header=T)
  
  # Filters
  synt_obs <- synt_obs[which(synt_obs$BLAT_match == 1 & synt_obs$align_score > 0.1),]
  synt_simul <- synt_simul[which(synt_simul$BLAT_match == 1 & synt_simul$align_score > 0.1),]
  synt_obs <- synt_obs %>% replace_with_na(replace = list(target_dist = "trans"))
  synt_simul <- synt_simul %>% replace_with_na(replace = list(target_dist = "trans"))
  
  if (sp == "opossum"){
    synt_obs <- synt_obs[which(!is.na(synt_obs$target_dist) & as.numeric(synt_obs$target_dist) < 2000000),]
    synt_simul <- synt_simul[which(!is.na(synt_simul$target_dist) & as.numeric(synt_simul$target_dist) < 2000000),]
    synt_obs_opposum <-  data.frame(ID = do.call(paste,c(synt_obs[c("origin_gene","origin_enh")],sep="-")))
    synt_simul_opposum <-  data.frame(ID = do.call(paste,c(synt_simul[c("origin_gene","origin_enh")],sep=":")))
    
  }else{
    # Distance classes
    max_dist = 2500000
    synt_obs$class_dist <-cut(synt_obs$origin_dist, breaks=seq(from=0, to=max_dist, by=50000), include.lowest = T)
    synt_simul$class_dist <-cut(synt_simul$origin_dist, breaks=seq(from=0, to=max_dist, by=50000), include.lowest = T)
    class_leg <- c("0", "500Kb", "1Mb", "1.5Mb", "2Mb", "2.5Mb")
    
    # Calculate proportion
    conserv <- data.frame(obs = sapply(levels(synt_obs$class_dist), function(x)
      nrow(synt_obs[which(!is.na(synt_obs$target_dist) & as.numeric(synt_obs$target_dist) < 2000000 & synt_obs$class_dist == x ),])
      /nrow(synt_obs[which(synt_obs$class_dist == x ),])))
    
    conserv$conf_low_obs <- sapply(levels(synt_obs$class_dist), function(x)
      prop.test(x = nrow(synt_obs[which(!is.na(synt_obs$target_dist) & as.numeric(synt_obs$target_dist) < 2000000 & synt_obs$class_dist == x ),]),
                n=nrow(synt_obs[which(synt_obs$class_dist == x ),]), p=0.5)$conf.int[1])
    
    conserv$conf_up_obs <-sapply(levels(synt_obs$class_dist), function(x)
      prop.test(x = nrow(synt_obs[which(!is.na(synt_obs$target_dist) & as.numeric(synt_obs$target_dist) < 2000000 & synt_obs$class_dist == x ),]),
                n=nrow(synt_obs[which(synt_obs$class_dist == x ),]), p=0.5)$conf.int[2])
    
    ### Intersection oppossum 
    synt_obs$ID <-  do.call(paste,c(synt_obs[c("origin_gene","origin_enh")],sep="-"))
    synt_simul$ID <- do.call(paste,c(synt_simul[c("origin_gene","origin_enh")],sep=":"))
    
    conserv$obs_opposum <- sapply(levels(synt_obs$class_dist), function(x)
      nrow(synt_obs[which(!is.na(synt_obs$target_dist) & as.numeric(synt_obs$target_dist) < 2000000 & synt_obs$class_dist == x & synt_obs$ID %in% synt_obs_opposum$ID),])
      /nrow(synt_obs[which(synt_obs$class_dist == x & synt_obs$ID %in% synt_obs_opposum$ID),]))
    
    conserv$conf_low_obs_opossum <- sapply(levels(synt_obs$class_dist), function(x)
      prop.test(x = nrow(synt_obs[which(!is.na(synt_obs$target_dist) & as.numeric(synt_obs$target_dist) < 2000000 & synt_obs$class_dist == x & synt_obs$ID %in% synt_obs_opposum$ID),]),
                n=nrow(synt_obs[which(synt_obs$class_dist == x & synt_obs$ID %in% synt_obs_opposum$ID),]), p=0.5)$conf.int[1])
    
    conserv$conf_up_obs_opossum  <-sapply(levels(synt_obs$class_dist), function(x)
      prop.test(x = nrow(synt_obs[which(!is.na(synt_obs$target_dist) & as.numeric(synt_obs$target_dist) < 2000000 & synt_obs$class_dist == x & synt_obs$ID %in% synt_obs_opposum$ID),]),
                n=nrow(synt_obs[which(synt_obs$class_dist == x & synt_obs$ID %in% synt_obs_opposum$ID),]), p=0.5)$conf.int[2])

    
  }

}

