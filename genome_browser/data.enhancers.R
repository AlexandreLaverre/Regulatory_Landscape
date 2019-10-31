#####################################################################

pathAnouk="/beegfs/data/necsulea/RegulatoryLandscapes_AN/"
pathAlex="/beegfs/data/alaverre/Regulatory_landscape/"

pathEnhancers=paste(pathAlex, "result/Supplementary_dataset3_enhancers/", sep="")

#####################################################################

enhancer.coords=list()

for(sp in c("human", "mouse")){
  print(sp)
  
  enhancer.coords[[sp]]=list()

  files=system(paste("ls ", pathEnhancers, sp, "/coord_enh/ | grep .bed",sep=""), intern=T)
  datasets=unlist(lapply(files, function(x) unlist(strsplit(x, split="_"))[1]))

  for(i in 1:length(files)){
    dataset=datasets[i]

    print(dataset)

    this.coords=read.table(paste(pathEnhancers, sp, "/coord_enh/", files[i], sep=""), h=T, stringsAsFactors=F, sep="\t")
    this.coords$chr=unlist(lapply(this.coords$chr, function(x) substr(x, 4, nchar(x))))
    this.coords=this.coords[, c("chr", "start", "end")]

    enhancer.coords[[sp]][[dataset]]=this.coords
  }
}

#####################################################################

save(enhancer.coords, file="RData/data.enhancers.RData")

#####################################################################

