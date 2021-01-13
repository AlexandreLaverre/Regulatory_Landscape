#######################################################################################
options(stringsAsFactors = FALSE)

objects=ls()

if(!"pathScripts"%in%objects){
  load=T
  
  source("parameters.R")
}

#######################################################################################

if(load){
  sp="human"
  load(paste(pathFigures, "RData/data.gene.enhancer.contacts.RData", sep=""))
  load(paste(pathFigures, "RData/data.enhancer.statistics.RData", sep="")) ## enhancers are already filtered for duplication levels, repeat proprtion etc
  load(paste(pathFigures, "RData/data.", sp, ".regland.conservation.RData",sep=""))
  load(paste(pathFigures, "RData/data.", sp, ".CM2019.SomaticOrgans.expdiv.Rdata", sep=""))
  load=FALSE
  
}

#######################################################################################
# Total contacts and median distance to enhancers

for (data in c("real", "simulated")){
  contact.data = gene.enhancer.contacts[["human"]][["ENCODE"]][[data]]
  
  # Filters on enhancers
  if (data == "real"){data.stats="original"}else{data.stats=data}
  enh.stat.data = enhancer.statistics[["human"]][["ENCODE"]][[data.stats]]
  contact.data = contact.data[which(contact.data$enhancer%in%enh.stat.data$enh),]
  
  if (data == "real"){data.name="original"}else{data.name=data}
  
  nb_total = unlist(with(contact.data, tapply(enhancer, factor(gene), function(x) length(x))))
  median_dist = unlist(with(contact.data, tapply(dist, factor(gene), median, na.rm=T)))
  
  gene.stat <- data.frame("nb_total"=nb_total, "median_dist"=median_dist)
  
  write.table(rownames(gene.stat[order(-gene.stat$nb_total),]), file =paste(pathFigures,"ranked.list.genes/ordered.genes.", data.name, ".total.contacts.txt", sep=""),
              row.names=FALSE, col.names=FALSE, quote=F)
  
  write.table(rownames(gene.stat[order(-gene.stat$median_dist),]), file =paste(pathFigures,"ranked.list.genes/ordered.genes.", data.name, ".median.enhancers.distance.txt", sep=""),
              row.names=FALSE, col.names=FALSE, quote=F)
}


#######################################################################################
# Regulatory Landscape Conservation 
ranked.var = c(".align.score.", ".conserved.seq.", ".conserved.synt.", ".conserv.contacts.")
names(ranked.var) = c("align_score", "ratio_cons_seq", "ratio_cons_synt", "ratio_cons_int")

for (data in c("obs", "sim")){
  regland = genes.conservation[["ENCODE"]][[data]][["all"]]
  
  if (data == "obs"){data.name="original"}else{data.name="simulated"}
  
  for (var in names(ranked.var)){
    write.table(rownames(regland[order(-regland[[var]]),]), file =paste(pathFigures,"ranked.list.genes/ordered.genes.", data.name, ranked.var[var], ".txt", sep=""), row.names=FALSE, col.names=FALSE, quote=F)
  }

}
###########

#######################################################################################
# Gene expression pattern Conservation 

write.table(rownames(expdiv[order(-expdiv$CorrelationSpearman),]), file =paste(pathFigures,"ranked.list.genes/ordered.genes.spearman.correlation.txt", sep=""), row.names=FALSE, col.names=FALSE, quote=F)
write.table(rownames(expdiv[order(-expdiv$CorrectedSpearman),]), file =paste(pathFigures,"ranked.list.genes/ordered.genes.corrected.spearman.correlation.txt", sep=""), row.names=FALSE, col.names=FALSE, quote=F)

write.table(rownames(expdiv[order(expdiv$EuclideanDistance),]), file =paste(pathFigures,"ranked.list.genes/ordered.genes.euclidean.distance.txt", sep=""), row.names=FALSE, col.names=FALSE, quote=F)
write.table(rownames(expdiv[order(-expdiv$CorrectedEuclideanSimilarity),]), file =paste(pathFigures,"ranked.list.genes/ordered.genes.corrected.euclidean.distance.txt", sep=""), row.names=FALSE, col.names=FALSE, quote=F)


###########
