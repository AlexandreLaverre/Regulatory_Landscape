###########################################################################

dirs=unlist(strsplit(getwd(), split="\\/"))
path=paste(dirs[1:(length(dirs)-2)], collapse="/")

pathFOCS=paste(path, "/data/FOCS/", sep="")

###########################################################################

sp="human"

load(paste(pathFOCS, sp, "/ENCODE/encode.enh.rpkm.RData", sep=""))
load(paste(pathFOCS, sp, "/ENCODE/encode.prom.rpkm.RData", sep=""))
version="hg19"


id.enhancers=rownames(Me_fpkm)
chr.enhancers=unlist(lapply(id.enhancers, function(x) unlist(strsplit(x, split=":"))[1]))
chr.enhancers=unlist(lapply(chr.enhancers, function(x) unlist(strsplit(x, split="_"))[2]))
pos.enhancers=unlist(lapply(id.enhancers, function(x) unlist(strsplit(x, split=":"))[2]))
start.enhancers=unlist(lapply(pos.enhancers, function(x) unlist(strsplit(x, split="_"))[1]))
end.enhancers=unlist(lapply(pos.enhancers, function(x) unlist(strsplit(x, split="_"))[2]))

results=data.frame("chr"=chr.enhancers, "start"=start.enhancers, "end"=end.enhancers, "id"=id.enhancers, stringsAsFactors=F)

write.table(results, file=paste(pathFOCS, sp, "/ENCODE/enhancer_coordinates_",version,".bed", sep=""), row.names=F, col.names=F, sep="\t", quote=F)

id.promoters=rownames(Mg_fpkm)
chr.promoters=unlist(lapply(id.promoters, function(x) unlist(strsplit(x, split=":"))[1]))
chr.promoters=unlist(lapply(chr.promoters, function(x) unlist(strsplit(x, split="_"))[2]))
pos.promoters=unlist(lapply(id.promoters, function(x) unlist(strsplit(x, split=":"))[2]))
start.promoters=unlist(lapply(pos.promoters, function(x) unlist(strsplit(x, split="_"))[1]))
end.promoters=unlist(lapply(pos.promoters, function(x) unlist(strsplit(x, split="_"))[2]))

results=data.frame("chr"=chr.promoters, "start"=start.promoters, "end"=end.promoters, "id"=id.promoters, stringsAsFactors=F)

write.table(results, file=paste(pathFOCS, sp, "/ENCODE/promoter_coordinates_",version,".bed", sep=""), row.names=F, col.names=F, sep="\t", quote=F)

  ## write expression data now

if(dim(Me_fpkm)[2]!=dim(Mg_fpkm)[2]){
  stop("weird! different column names")
}

if(!all(colnames(Me_fpkm)%in%colnames(Mg_fpkm))){
  stop("weird! different column names")
}

Me_fpkm=Me_fpkm[,colnames(Mg_fpkm)]

write.table(Me_fpkm, file=paste(pathFOCS, sp, "/ENCODE/enhancer_activity.txt", sep=""), row.names=T, col.names=T, sep="\t", quote=F)

write.table(Mg_fpkm, file=paste(pathFOCS, sp, "/ENCODE/promoter_activity.txt", sep=""), row.names=T, col.names=T, sep="\t", quote=F)


###########################################################################


