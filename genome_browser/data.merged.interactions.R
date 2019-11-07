#####################################################################

pathAnouk="/beegfs/data/necsulea/RegulatoryLandscapes_AN/"
pathAlex="/beegfs/data/alaverre/Regulatory_landscape/"

pathInteractions=paste(pathAlex, "result/Supplementary_dataset2_merged_fragments/", sep="")

#####################################################################

maxdist=10e6 ## 10Mb maximum distance between bait and contacted position

#####################################################################

merged.interactions=list()

for(sp in c("human", "mouse")){

  this.mi=read.table(paste(pathInteractions, sp,"/all_interactions_merged.txt", sep=""), h=T, stringsAsFactors=F, sep="\t")
  this.mi$bait_chr=unlist(lapply(this.mi$bait_chr, function(x) substr(x, 4, nchar(x))))
  this.mi$chr=unlist(lapply(this.mi$chr, function(x) substr(x, 4, nchar(x))))

  this.mi$bait_midpos=(this.mi$bait_start+this.mi$bait_end)/2
  this.mi$otherEnd_midpos=(this.mi$start+this.mi$end)/2
  this.mi$distance=abs(this.mi$bait_midpos-this.mi$otherEnd_midpos)

  this.mi=this.mi[which(this.mi$distance<=maxdist),]
  
  this.mi$bait_ID=paste(this.mi$bait_chr, this.mi$bait_start, this.mi$bait_end, sep=",")

  if(all(this.mi$bait_chr==this.mi$chr)){
    print("all in cis")
  }

  merged.interactions[[sp]]=this.mi
}

#####################################################################

save(merged.interactions, file="RData/data.merged.interactions.RData")

#####################################################################
