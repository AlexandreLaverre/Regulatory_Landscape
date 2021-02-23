#######################################################################################

options(stringsAsFactors = FALSE)

#######################################################################################

objects=ls()

if(!"pathScripts"%in%objects){
  load=T
  
  enh="ENCODE"
  
  source("../main_figures/parameters.R")
}

#######################################################################################

for(sp in c("human", "mouse")){

  if(sp=="human"){
    spname="Human"
  }

  if(sp=="mouse"){
    spname="Mouse"
  }
 
  load(paste(pathFigures, "RData/data.gene.enhancer.contacts.RData", sep=""))
  contact.data = gene.enhancer.contacts[[sp]][[enh]][["real"]]
    
  load(paste(pathFigures, "RData/data.regland.conservation.RData",sep=""))
  regcons=regland.conservation[[sp]][[enh]]

  load(paste(pathFigures, "RData/data.", sp, ".CM2019.SomaticOrgans.expdiv.RData", sep=""))
 
  pathResults=paste(pathFinalData, "SupplementaryDataset9/", sep="")

  
  ## Total contacts and median distance to enhancers
  ## only PCHi-C data
  
  ## total number of contacted enhancers
  nb.total=tapply(contact.data$enhancer, as.factor(contact.data$gene), function(x) length(unique(x)))
  
  df.nbcontacts=data.frame("GeneID"=levels(as.factor(contact.data$gene)), "NbContacts"=nb.total, stringsAsFactors=F)
  df.nbcontacts=df.nbcontacts[order(df.nbcontacts$NbContacts, decreasing=T),]
  
  write.table(df.nbcontacts, file=paste(pathResults, "GeneTable_TotalNbContacts_", spname,".txt", sep=""), row.names=F, col.names=T, sep="\t", quote=F)
  writeLines(df.nbcontacts$GeneID, con=paste(pathResults, "GeneList_TotalNbContacts_", spname,".txt", sep=""))
  
  ## median distance
  median.dist = tapply(contact.data$dist, as.factor(contact.data$gene), median, na.rm=T)
  
  df.mediandist=data.frame("GeneID"=levels(as.factor(contact.data$gene)), "MedianDistance"=median.dist, stringsAsFactors=F)
  df.mediandist=df.mediandist[order(df.mediandist$MedianDistance, decreasing=T),]
  
  write.table(df.mediandist, file=paste(pathResults, "GeneTable_MedianDistance_", spname,".txt", sep=""), row.names=F, col.names=T, sep="\t", quote=F)
  writeLines(df.mediandist$GeneID, con=paste(pathResults, "GeneList_MedianDistance_", spname,".txt", sep=""))
  
#######################################################################################
  
  ## enhancer alignment score
  df.alnscore=regcons[,c("gene", "mean.aln.score.all")]
  colnames(df.alnscore)=c("GeneID", "MeanAlignmentScore")
  
  df.alnscore=df.alnscore[which(!is.na(df.alnscore$MeanAlignmentScore)),]
  
  df.alnscore=df.alnscore[order(df.alnscore$MeanAlignmentScore, decreasing=T),]
  
  write.table(df.alnscore, file=paste(pathResults, "GeneTable_MeanEnhancerAlignmentScore_", spname,".txt", sep=""), row.names=F, col.names=T, sep="\t", quote=F)
  writeLines(df.alnscore$GeneID, con=paste(pathResults, "GeneList_MeanEnhancerAlignmentScore_", spname,".txt", sep=""))
  
  #######################################################################################
  ## contact conservation
  df.contactcons=regcons[,c("gene", "fr.contact.cons.all")]
  colnames(df.contactcons)=c("GeneID", "FractionConservedContacts")
  df.contactcons=df.contactcons[which(!is.na(df.contactcons$FractionConservedContacts)),]
  
  df.contactcons=df.contactcons[order(df.contactcons$FractionConservedContacts, decreasing=T),]
  
  write.table(df.contactcons, file=paste(pathResults, "GeneTable_FractionConservedContacts_", spname,".txt", sep=""), row.names=F, col.names=T, sep="\t", quote=F)
  writeLines(df.contactcons$GeneID, con=paste(pathResults, "GeneList_FractionConservedContacts_", spname,".txt", sep=""))
  
 #######################################################################################
  
  ## synteny conservation
  df.syntenycons=regcons[,c("gene", "fr.synteny.cons.all")]
  colnames(df.syntenycons)=c("GeneID", "FractionConservedSynteny")
  df.syntenycons=df.syntenycons[which(!is.na(df.syntenycons$FractionConservedSynteny)),]
  
  df.syntenycons=df.syntenycons[order(df.syntenycons$FractionConservedSynteny, decreasing=T),]
  
  write.table(df.syntenycons, file=paste(pathResults, "GeneTable_FractionConservedSynteny_", spname,".txt", sep=""), row.names=F, col.names=T, sep="\t", quote=F)
  writeLines(df.syntenycons$GeneID, con=paste(pathResults, "GeneList_FractionConservedSynteny_", spname,".txt", sep=""))
  
 #######################################################################################
  
  ## expression conservation
  df.expcons=expdiv[,c("IDHuman", "CorrectedSpearman")]
  colnames(df.expcons)=c("GeneID", "CorrectedSpearmanCorrelation")
  
  df.expcons=df.expcons[which(!is.na(df.expcons$CorrectedSpearmanCorrelation)),]
  df.expcons=df.expcons[order(df.expcons$CorrectedSpearmanCorrelation, decreasing=T),]
  
  write.table(df.expcons, file=paste(pathResults, "GeneTable_ExpressionConservation_", spname,".txt", sep=""), row.names=F, col.names=T, sep="\t", quote=F)
  writeLines(df.expcons$GeneID, con=paste(pathResults, "GeneList_ExpressionConservation_", spname,".txt", sep=""))
  
 #######################################################################################

}

#######################################################################################
