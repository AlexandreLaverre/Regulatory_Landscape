###########################################################################

dirs=unlist(strsplit(getwd(), split="\\/"))
path=paste(dirs[1:(length(dirs)-2)], collapse="/")

pathFOCS=paste(path, "/data/FOCS/", sep="")

###########################################################################

sp="human"

load(paste(pathFOCS, sp, "/GRO-seq/groseq.enh.rpkm.RData", sep=""))
load(paste(pathFOCS, sp, "/GRO-seq/groseq.prom.rpkm.RData", sep=""))
load(paste(pathFOCS, sp, "/GRO-seq/groseq.enh.pos.RData", sep=""))
load(paste(pathFOCS, sp, "/GRO-seq/groseq.prom.pos.RData", sep=""))

version="hg19"

enh.bs=as.data.frame(enh.bs)
prom.bs=as.data.frame(prom.bs)

print(all(prom.bs$ent_id==rownames(Mg_rpkm)))

## the rownames on Me_rpkm are not informative! I will assume that they are in the same order as the rows in enh.bs

rownames(Me_rpkm)=enh.bs$name

id.enhancers=rownames(Me_rpkm)
chr.enhancers=unlist(lapply(id.enhancers, function(x) unlist(strsplit(x, split=":"))[1]))
chr.enhancers=unlist(lapply(chr.enhancers, function(x) unlist(strsplit(x, split="_"))[2]))
pos.enhancers=unlist(lapply(id.enhancers, function(x) unlist(strsplit(x, split=":"))[2]))
start.enhancers=unlist(lapply(pos.enhancers, function(x) unlist(strsplit(x, split="_"))[1]))
end.enhancers=unlist(lapply(pos.enhancers, function(x) unlist(strsplit(x, split="_"))[2]))

results=data.frame("chr"=chr.enhancers, "start"=start.enhancers, "end"=end.enhancers, "id"=id.enhancers, stringsAsFactors=F)

write.table(results, file=paste(pathFOCS, sp, "/GRO-seq/enhancer_coordinates_",version,".bed", sep=""), row.names=F, col.names=F, sep="\t", quote=F)


## promoter coordinates


id.promoters=rownames(prom.bs)

chr.promoters=prom.bs$seqnames
strand.promoters=prom.bs$strand

start.promoters=rep(NA, dim(prom.bs)[1])
end.promoters=rep(NA, dim(prom.bs)[1])

promsize=1000
start.promoters[which(prom.bs$strand=="+")]=prom.bs$start[which(prom.bs$strand=="+")]-promsize
end.promoters[which(prom.bs$strand=="+")]=prom.bs$start[which(prom.bs$strand=="+")]+promsize

start.promoters[which(prom.bs$strand=="-")]=prom.bs$end[which(prom.bs$strand=="-")]-promsize
end.promoters[which(prom.bs$strand=="-")]=prom.bs$end[which(prom.bs$strand=="-")]+promsize

results=data.frame("chr"=chr.promoters, "start"=start.promoters, "end"=end.promoters, "id"=id.promoters, "strand"=strand.promoters, stringsAsFactors=F)

write.table(results, file=paste(pathFOCS, sp, "/GRO-seq/promoter_coordinates_",version,".bed", sep=""), row.names=F, col.names=F, sep="\t", quote=F)

## write expression data now

if(dim(Me_rpkm)[2]!=dim(Mg_rpkm)[2]){
  stop("weird! different column names")
}

if(!all(colnames(Me_rpkm)%in%colnames(Mg_rpkm))){
  stop("weird! different column names")
}

Me_rpkm=Me_rpkm[,colnames(Mg_rpkm)]

write.table(Me_rpkm, file=paste(pathFOCS, sp, "/GRO-seq/enhancer_activity.txt", sep=""), row.names=T, col.names=T, sep="\t", quote=F)

write.table(Mg_rpkm, file=paste(pathFOCS, sp, "/GRO-seq/promoter_activity.txt", sep=""), row.names=T, col.names=T, sep="\t", quote=F)

###########################################################################


