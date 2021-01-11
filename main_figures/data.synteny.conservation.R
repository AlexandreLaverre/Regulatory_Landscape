####################################################################################################

library(data.table)

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
      
      synt_obs <- fread(paste(path_evol,"/synteny_conservation/", enh, "/", ref_sp, "2", sp, "_original_synteny.txt", sep=""), header=T)
      synt_simul <- fread(paste(path_evol,"/synteny_conservation/", enh, "/", ref_sp, "2", sp, "_simulated_synteny.txt", sep=""), header=T)

      class(synt_obs)<-"data.frame"
      class(synt_simul)<-"data.frame"
      
      ## sequences have to be aligned
      synt_obs <- synt_obs[which(synt_obs$align_score>0),]
      synt_simul <- synt_simul[which(synt_simul$align_score>0),]
      
      ## threshold alignment score: 5% quantile, observed values

      align.threshold=quantile(synt_obs$align, p=0.1)

      print(paste("alignment score threshold", align.threshold))
      
      ## enhancer have to be unduplicated and aligned at least at align.threshold in target genome
      synt_obs <- synt_obs[which(synt_obs$BLAT_match == 1 & synt_obs$align_score >= align.threshold),]
      synt_simul <- synt_simul[which(synt_simul$BLAT_match == 1 & synt_simul$align_score >=align.threshold),]

      ## compute distance classes for original pairs

      synt_obs$class_dist <- cut(synt_obs$origin_dist, breaks=seq(from=minDistance, to=maxDistance, by=50000), include.lowest = T)
      synt_simul$class_dist <- cut(synt_simul$origin_dist, breaks=seq(from=minDistance, to=maxDistance, by=50000), include.lowest = T)
      
      ## if in trans, target distance is NA
      synt_obs$target_dist[which(synt_obs$target_dist=="trans")] <- NA
      synt_simul$target_dist[which(synt_simul$target_dist=="trans")] <- NA 

      ## transform distance to numeric values
      synt_obs$target_dist <- as.numeric(synt_obs$target_dist)
      synt_simul$target_dist <- as.numeric(synt_simul$target_dist)
      
      conserv_synteny[[enh]][[sp]]=list("synt_obs"=synt_obs, "synt_simul"=synt_simul)
            
      message(sp, " done !")
    }
  }

  ## save data 
  
  #save(conserv_synteny, file = paste(pathFigures, "/RData/data.synteny.conservation.", ref_sp, ".RData", sep=""))
  
  ## all done!
}

#######################################################################################################################################
