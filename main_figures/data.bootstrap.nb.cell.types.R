###########################################################################################

source("parameters.R") ## paths are defined based on the user name

set.seed(19)

library(bootBCa, lib=pathRlibs)

###########################################################################################

for(sp in c("human", "mouse")){

  load(paste(pathFigures, "RData/data.fragment.contacts.RData", sep=""))
  load(paste(pathFigures, "RData/data.sample.info.RData", sep=""))

   obs=observed.contacts[[sp]]
  sim=simulated.contacts[[sp]]

  info=sampleinfo[[sp]]
  rownames(info)=info$Sample.ID
  
  samples=info$Sample.ID 
  celltypes=info$Broad.cell.type.or.tissue
  names(celltypes)=samples
  
  print(paste("there are", length(samples), "samples"))
  
  obs$dist_class <- cut(obs$distance, breaks=seq(from=minDistance, to=maxDistance, by=50e3), include.lowest = T)

  obs$nb_celltypes <- apply(obs[,samples],1, function(x) length(unique(celltypes[which(!is.na(x))])))
  obs$celltype_class<- cut(obs$nb_celltypes, breaks=c(0:7, max(obs$nb_celltypes)), include.lowest=T)
  levels(obs$celltype_class)=c(as.character(1:7), ">7")
  
  sim$dist_class <- cut(sim$distance, breaks=seq(from=minDistance, to=maxDistance, by=50e3), include.lowest = T)

  sim$nb_celltypes <- apply(sim[,samples],1, function(x) length(unique(celltypes[which(!is.na(x))])))
  sim$celltype_class<- cut(sim$nb_celltypes, breaks=c(0:7, max(sim$nb_celltypes)), include.lowest=T)
  levels(sim$celltype_class)=c(as.character(1:7), ">7")

  filtered_data <- list("Original"=obs, "Simulated"=sim) ## data is already unbaited, in cis, in the right distance range
 
  ## compute number of interactions in each nb samples class
  
  ## nb cell types
  
  nb_celltypes_matrix <- sapply(filtered_data, function(x) as.numeric(table(x$celltype_class)))
  rownames(nb_celltypes_matrix) = levels(obs$celltype_class)

  nb_celltypes_matrix=t(nb_celltypes_matrix)
  pc_nb_celltypes_matrix=100*nb_celltypes_matrix/apply(nb_celltypes_matrix,1, sum)

  ## divide interactions based on distance, compute mean number of samples by distance class

  print("computing bootstrap confidence intervals")

  dist_cons_celltypes <- lapply(filtered_data, function(x) tapply(x$nb_celltypes, as.factor(x$dist_class), function(y) {z<-BCa(y, delta=NA, M=100, theta=mean); return(z)}))
  
  dist_conf_low_celltypes <- t(sapply(dist_cons_celltypes, function(x) unlist(lapply(x, function(y) y[4]))))
  dist_conf_high_celltypes <- t(sapply(dist_cons_celltypes, function(x) unlist(lapply(x, function(y) y[5]))))

  mean_nb_celltypes_dist <- t(sapply(dist_cons_celltypes, function(x)   unlist(lapply(x, function(y) y[3])))) 

  mean_dist <- t(sapply(filtered_data, function(x)   tapply(x$distance, as.factor(x$dist_class), mean, na.rm=T)))

  print("done")

  ## save objects

  save(list=ls(), file=paste(pathFigures, "RData/data.bootstrap.nb.cell.types.",sp,".RData",sep=""))
  
}

############################################################################################
