####################################################################################################

options(stringsAsFactors = FALSE)

source("parameters.R")

all.species <- c("human", "macaque", "dog", "cow", "elephant", "rabbit", "mouse", "rat", "opossum", "chicken")

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
    conserv_synteny[[enh]] <- list()

    message("Running conserv synteny with ", enh)
    
    for (sp in species){
      conserv_synteny[[enh]][[sp]]<-list()
      
      synt_obs <- read.table(paste(path_evol,"/synteny_conservation/", enh, "/", ref_sp, "2", sp, "_original_synteny.txt", sep=""), header=T)
      synt_simul <- read.table(paste(path_evol,"/synteny_conservation/", enh, "/", ref_sp, "2", sp, "_simulated_synteny.txt", sep=""), header=T)
      
      ## enhancer filters
      
      synt_obs <- synt_obs[which(synt_obs$BLAT_match < 2 & synt_obs$align_score > 0.4 & synt_obs$origin_dist > minDistance & synt_obs$origin_dist < maxDistance),]
      synt_simul <- synt_simul[which(synt_simul$BLAT_match < 2 & synt_simul$align_score > 0.4 & synt_simul$origin_dist > minDistance & synt_simul$origin_dist < maxDistance),]

      ## compute distance classes for original pairs

      synt_obs$class_dist <- cut(synt_obs$origin_dist, breaks=seq(from=minDistance, to=maxDistance, by=50000), include.lowest = T)
      synt_simul$class_dist <- cut(synt_simul$origin_dist, breaks=seq(from=minDistance, to=maxDistance, by=50000), include.lowest = T)
      
      ## if in trans, target distance is NA
      synt_obs$target_dist[which(synt_obs$target_dist=="trans")] <- NA
      synt_simul$target_dist[which(synt_simul$target_dist=="trans")] <- NA 

      ## transform distance to numeric values
      synt_obs$target_dist <- as.numeric(synt_obs$target_dist)
      synt_simul$target_dist <- as.numeric(synt_simul$target_dist)
      
      ## calculate proportion of pairs conserved in synteny
      nb_cons_synt_obs <- length(which(synt_obs$target_dist < maxDistance))
      nb_cons_synt_simul <- length(which(synt_simul$target_dist < maxDistance))
      
      mat <- matrix(c(nb_cons_synt_obs, nb_cons_synt_simul, nrow(synt_obs)-nb_cons_synt_obs, nrow(synt_simul)-nb_cons_synt_simul), nrow=2)
      
      prop.test.obs <- prop.test(x = nb_cons_synt_obs, n=nrow(synt_obs), p=0.5)
      prop.test.simul <- prop.test(x = nb_cons_synt_simul, n=nrow(synt_simul), p=0.5)

      conserv_synteny[[enh]][[sp]]=list("synt_obs"=synt_obs, "synt_simul"=synt_simul, "prop_cons_obs"=(nb_cons_synt_obs/nrow(synt_obs)), "prop_cons_simul"=(nb_cons_synt_simul/nrow(synt_simul)), "conf_int_obs"=prop.test.obs$conf.int, "conf_int_simul"=prop.test.obs$conf.int, "pval"=prop.test(mat)$p.value)
      
      
      message(sp, " done !")
    }
  }

  ## save data 
  
  save(conserv_synteny, file = paste(pathFigures, "/RData/data.synteny.conservation.", ref_sp, ".RData", sep=""))
  
  ## all done!
}

#######################################################################################################################################
