#############################################################

pathRoadmap="../../data/RoadmapEpigenomics/hg19/regulatory_regions/"

#############################################################

load(paste(pathRoadmap,"promoter_regions/promoter.positions.FOCS.RData", sep=""))
load(paste(pathRoadmap,"enhancer_regions/enhancer.positions.FOCS.RData", sep=""))

#############################################################

prom=as.data.frame(prom.bs)

prom=prom[,c("seqnames","start","end")]
prom$id=paste(prom$seqnames, prom$start, prom$end, sep=",")

write.table(prom, file=paste(pathRoadmap,"promoter_regions/promoters_FOCS.bed",sep=""), row.names=F, col.names=F,sep="\t", quote=F)

#############################################################

enh=as.data.frame(enh.bs)

enh=enh[,c("seqnames","start","end")]
enh$id=paste(enh$seqnames, enh$start, enh$end, sep=",")

write.table(enh, file=paste(pathRoadmap,"enhancer_regions/enhancers_FOCS.bed",sep=""), row.names=F, col.names=F,sep="\t", quote=F)

#############################################################
