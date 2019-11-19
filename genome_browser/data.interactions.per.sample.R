#####################################################################

pathAnouk="/beegfs/data/necsulea/RegulatoryLandscapes_AN/"
pathAlex="/beegfs/data/alaverre/Regulatory_landscape/"

pathInteractions=paste(pathAlex, "result/Supplementary_dataset1_original_fragments/", sep="")

#####################################################################

interactions=list()

for(sp in c("human", "mouse")){
  
  interactions[[sp]]=list()

  files=system(paste("ls ", pathInteractions, sp, "/interactions_samples/ | grep ibed", sep=""), intern=T)
  
  samples=unlist(lapply(files, function(x) unlist(strsplit(x, split="\\."))[1]))
  
  for(sample in samples){
    data=read.table(paste(pathInteractions, sp, "/interactions_samples/", sample,".ibed", sep=""), h=T, stringsAsFactors=F)

    ## select cis interactions

    data=data[which(data$bait_chr==data$otherEnd_chr),]

    interactions[[sp]][[sample]]=data
   
  }
}

#####################################################################

save(interactions, file="RData/data.interactions.per.sample.RData")

#####################################################################
