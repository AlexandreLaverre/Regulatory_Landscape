#####################################################################

pathAnouk="/beegfs/data/necsulea/RegulatoryLandscapes_AN/"
pathAlex="/beegfs/data/alaverre/Regulatory_landscape/"

pathInteractions=paste(pathAlex, "result/Supplementary_dataset2_merged_fragments/", sep="")

#####################################################################

merged.interactions=list()

for(sp in c("human", "mouse")){

  this.mi=read.table(paste(pathInteractions, sp,"/all_interactions_merged.txt", sep=""), h=T, stringsAsFactors=F, sep="\t")
  this.mi$bait_chr=unlist(lapply(this.mi$bait_chr, function(x) substr(x, 3, nchar(x))))
  this.mi$chr=unlist(lapply(this.mi$chr, function(x) substr(x, 3, nchar(x))))
  
  this.mi$bait_ID=paste(this.mi$bait_chr, this.mi$bait_start, this.mi$bait_end, sep=",")

  if(all(this.mi$bait_chr==this.mi$chr)){
    print("all in cis")
  }

  merged.interactions[[sp]]=this.mi
}

#####################################################################

save(merged.interactions, file="RData/data.merged.interaction.RData")

#####################################################################
