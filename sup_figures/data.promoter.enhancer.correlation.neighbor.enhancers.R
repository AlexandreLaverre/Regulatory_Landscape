#########################################################################################################################

source("../main_figures/parameters.R") ##  are defined based on the user name

library(Hmisc)
library(data.table)
library(bootBCa, lib=pathRlibs)

set.seed(19)

options(stringsAsFactors = FALSE)

load(paste(pathFigures,  "RData/data.fragment.contacts.RData",sep=""))

thisMaxDistance=1e6

##################################  Correlation gene expression & enhancers activity ####################################

for(ref_sp in c("human", "mouse")){

  ## fragment contacts
  
  frag.contact.obs=observed.contacts[[ref_sp]]
  
  frag.contact.obs$idcontact=paste(frag.contact.obs$id_bait, frag.contact.obs$id_frag, sep="-")
    
  enhancers = enhancer.datasets[[ref_sp]]
  
  obs_correl_activity_dist <- list()
  neighbors_correl_activity_dist <- list()
    
  for (enh in enhancers){
    if(file.exists(paste(pathFinalData, "SupplementaryDataset8", ref_sp, enh, "expression_correlations_real_data.txt", sep="/")) & file.exists(paste(pathFinalData, "SupplementaryDataset8", ref_sp, enh, "expression_correlations_neighbor_promoters_enhancers.txt", sep="/"))){

      print(paste(ref_sp, enh))
      
      obs <- fread(paste(pathFinalData, "SupplementaryDataset8", ref_sp, enh, "expression_correlations_real_data.txt", sep="/"), h=T, sep="\t")
      class(obs)<-"data.frame"

      ## keep only previously filtered bait-fragment contacts
      obs$IDBait=unlist(lapply(obs$IDBait, function(x) paste(unlist(strsplit(x, split=",")), collapse=":")))
      obs$IDContactedFragment=unlist(lapply(obs$IDContactedFragment, function(x) paste(unlist(strsplit(x, split=",")),collapse=":")))
      obs$IDContact=paste(obs$IDBait, obs$IDContactedFragment, sep="-")

      obs=obs[which(obs$IDContact%in%frag.contact.obs$idcontact),]
      
      neighbors <- fread(paste(pathFinalData, "SupplementaryDataset8", ref_sp, enh, "expression_correlations_neighbor_promoters_enhancers.txt", sep="/"), h=T, sep="\t")
      class(neighbors)<-"data.frame"
            
      # distance classes
    
      obs$dist_class <-cut(obs$Distance, breaks=seq(from=minDistance, to=thisMaxDistance+25000, by=25000), include.lowest = T)
      neighbors$dist_class <- cut(neighbors$Distance, breaks=seq(from=minDistance, to=thisMaxDistance+25000, by=25000), include.lowest = T)

      ## same promoters in each distance class
      common.class=lapply(levels(obs$dist_class), function(x) intersect(obs$IDPromoter[which(obs$dist_class==x)], neighbors$IDPromoter[which(neighbors$dist_class==x)]))
      names(common.class)=levels(obs$dist_class)

      nb.prom=unlist(lapply(common.class, length))
      common.class=common.class[which(nb.prom>=10)]

      levels=levels(obs$dist_class)
      levels=levels[which(levels%in%names(common.class))]

      BC.obs=lapply(levels, function(x) {y=obs$SpearmanCorrelation[which(obs$dist_class==x & obs$IDPromoter%in%common.class[[x]])]; return(BCa(y, delta=NA, M=100, theta=mean, na.rm=T))})
    
      obs_correl_activity_dist[[enh]]=unlist(lapply(BC.obs, function(x) x[3]))
      names(obs_correl_activity_dist[[enh]])=levels
      obs_correl_activity_dist[[paste0(enh, "_conflow")]]=unlist(lapply(BC.obs, function(x) x[4]))
      obs_correl_activity_dist[[paste0(enh, "_confup")]]=unlist(lapply(BC.obs, function(x) x[5]))
       
      BC.neighbors=lapply(levels, function(x) {y=neighbors$SpearmanCorrelation[which(neighbors$dist_class==x & neighbors$IDPromoter%in%common.class[[x]])]; return(BCa(y, delta=NA, M=100, theta=mean, na.rm=T))})
      neighbors_correl_activity_dist[[enh]]=unlist(lapply(BC.neighbors, function(x) x[3]))
      names(obs_correl_activity_dist[[enh]])=levels
      neighbors_correl_activity_dist[[paste0(enh, "_conflow")]]=unlist(lapply(BC.neighbors, function(x) x[4]))
      neighbors_correl_activity_dist[[paste0(enh, "_confup")]]=unlist(lapply(BC.neighbors, function(x) x[5]))

    }
  }
  
  correl_activity <- list(obs=obs_correl_activity_dist, neighbors=neighbors_correl_activity_dist)
  
################################################# Save RData ################################################# 
  
  save(correl_activity, file = paste(pathFigures, "RData/data.promoter.enhancer.correlation.neighbor.enhancers.", ref_sp, ".RData", sep=""))
 
 #########################################################################################################################
}

#########################################################################################################################
