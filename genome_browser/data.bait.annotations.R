#####################################################################

pathAnouk="/beegfs/data/necsulea/RegulatoryLandscapes_AN/"
pathAlex="/beegfs/data/alaverre/Regulatory_landscape/"

pathInteractions=paste(pathAlex, "result/Supplementary_dataset1_original_fragments/", sep="")

#####################################################################

bait.overlapTSS=list()

for(sp in c("human", "mouse")){ 
 
  ovtss=read.table(paste(pathInteractions, sp, "/bait_overlap_TSS_1kb.txt", sep=""), h=T, stringsAsFactors=F, sep="\t")

  ovtss$chr=unlist(lapply(ovtss$chr, function(x) substr(x, 4, nchar(x))))
  ovtss$ID=paste(ovtss$chr, ovtss$start, ovtss$end, sep=",")
  
  ovtss=ovtss[which(!is.na(ovtss$gene_ID)),]
  
  ovtss$nbgenes=unlist(lapply(ovtss$gene_ID, function(x) length(unlist(strsplit(x, split=",")))))
  

  all.bait.ids=rep(ovtss$fragment_ID, ovtss$nbgenes)
  all.gene.ids=unlist(lapply(ovtss$gene_ID, function(x) unlist(strsplit(x, split=","))))

  stop()

  bait.overlapTSS[[sp]]=data.frame("gene_ID"=all.gene.ids, "bait_ID"=all.bait.ids, stringsAsFactors=F)

}

#####################################################################

save(annot.baits.TSS, file="RData/data.interactions.annotations.RData")

#####################################################################
