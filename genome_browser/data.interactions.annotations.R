#####################################################################

pathAnouk="/beegfs/data/necsulea/RegulatoryLandscapes/"
pathAlex="/beegfs/data/alaverre/Regulatory_landscape/"

pathInteractions=paste(pathAlex, "result/Supplementary_dataset1_original_fragments/", sep="")

#####################################################################

annot.baits.TSS=list()

for(sp in c("human", "mouse")){ 
  bait.coords=read.table(paste(pathInteractions, sp, "/bait_coord.txt", sep=""), h=T, stringsAsFactors=F, sep="\t")
  bait.coords$fragment_ID=paste(bait.coords$chr, bait.coords$start, bait.coords$end, sep=",")
    
  frag.overlapTSS=read.table(paste(pathInteractions, sp, "/fragments_overlap_TSS.txt", sep=""), h=T, stringsAsFactors=F, sep="\t")
  frag.overlapTSS$fragment_ID=paste(frag.overlapTSS$chr, frag.overlapTSS$start, frag.overlapTSS$end, sep=",")
  bait.overlapTSS=frag.overlapTSS[which(frag.overlapTSS$fragment_ID%in%bait.coords$fragment_ID),]
  
  bait.overlapTSS=bait.overlapTSS[which(!is.na(bait.overlapTSS$gene_ID)),]
  bait.overlapTSS$nbgenes=unlist(lapply(bait.overlapTSS$gene_ID, function(x) length(unlist(strsplit(x, split=",")))))
  

  all.bait.ids=rep(bait.overlapTSS$fragment_ID, bait.overlapTSS$nbgenes)
  all.gene.ids=unlist(lapply(bait.overlapTSS$gene_ID, function(x) unlist(strsplit(x, split=","))))

  annot.baits.TSS[[sp]]=data.frame("gene_ID"=all.gene.ids, "bait_ID"=all.bait.ids, stringsAsFactors=F)

}

#####################################################################

save(annot.baits.TSS, file="RData/data.interactions.annotations.RData")

#####################################################################
