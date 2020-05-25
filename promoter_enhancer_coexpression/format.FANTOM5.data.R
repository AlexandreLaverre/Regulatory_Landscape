###########################################################################

dirs=unlist(strsplit(getwd(), split="\\/"))
path=paste(dirs[1:(length(dirs)-2)], collapse="/")

pathFOCS=paste(path, "/data/FOCS/", sep="")

###########################################################################

for(sp in c("human", "mouse")){
  if(sp=="human"){
    load(paste(pathFOCS, sp, "/FANTOM5/fantom.enh.tpm.RData", sep=""))
    load(paste(pathFOCS, sp, "/FANTOM5/fantom.prom.tpm.RData", sep=""))
    version="hg19"
    
  }

  if(sp=="mouse"){
    load(paste(pathFOCS, sp, "/FANTOM5/fantom.enh.tpm.mm9.RData", sep=""))
    load(paste(pathFOCS, sp, "/FANTOM5/fantom.prom.tpm.mm9.RData", sep=""))
    version="mm9"
  }

  id.enhancers=rownames(Me_fpkm)
  chr.enhancers=unlist(lapply(id.enhancers, function(x) unlist(strsplit(x, split=":"))[1]))
  pos.enhancers=unlist(lapply(id.enhancers, function(x) unlist(strsplit(x, split=":"))[2]))
  start.enhancers=unlist(lapply(pos.enhancers, function(x) unlist(strsplit(x, split="-"))[1]))
  end.enhancers=unlist(lapply(pos.enhancers, function(x) unlist(strsplit(x, split="-"))[2]))

  results=data.frame("chr"=chr.enhancers, "start"=start.enhancers, "end"=end.enhancers, "id"=id.enhancers, stringsAsFactors=F)

  write.table(results, file=paste(pathFOCS, sp, "/FANTOM5/enhancer_coordinates_",version,".bed", sep=""), row.names=F, col.names=F, sep="\t", quote=F)

  id.promoters=rownames(Mg_fpkm)
  chr.promoters=unlist(lapply(id.promoters, function(x) unlist(strsplit(x, split=":"))[1]))
  pos.promoters=unlist(lapply(id.promoters, function(x) unlist(strsplit(x, split=":"))[2]))
  strand.promoters=unlist(lapply(pos.promoters, function(x) unlist(strsplit(x, split=","))[2]))
  pos.promoters=unlist(lapply(pos.promoters, function(x) unlist(strsplit(x, split=","))[1]))
  start.promoters=unlist(lapply(pos.promoters, function(x) unlist(strsplit(x, split="\\.\\."))[1]))
  end.promoters=unlist(lapply(pos.promoters, function(x) unlist(strsplit(x, split="\\.\\."))[2]))

  results=data.frame("chr"=chr.promoters, "start"=start.promoters, "end"=end.promoters, "id"=id.promoters, "strand"=strand.promoters, stringsAsFactors=F)

  write.table(results, file=paste(pathFOCS, sp, "/FANTOM5/promoter_coordinates_",version,".bed", sep=""), row.names=F, col.names=F, sep="\t", quote=F)

  ## write expression data now

  if(dim(Me_fpkm)[2]!=dim(Mg_fpkm)[2]){
    stop("weird! different column names")
  }

  if(!all(colnames(Me_fpkm)%in%colnames(Mg_fpkm))){
    stop("weird! different column names")
  }

  Me_fpkm=Me_fpkm[,colnames(Mg_fpkm)]
  
  write.table(Me_fpkm, file=paste(pathFOCS, sp, "/FANTOM5/enhancer_activity.txt", sep=""), row.names=T, col.names=T, sep="\t", quote=F)
  
  write.table(Mg_fpkm, file=paste(pathFOCS, sp, "/FANTOM5/promoter_activity.txt", sep=""), row.names=T, col.names=T, sep="\t", quote=F)
  
   
  ## cleanup
  rm("Me_fpkm")
  rm("Mg_fpkm")
}

###########################################################################


