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

  if(length(files)>0){
    datasets=unlist(lapply(files, function(x) unlist(strsplit(x, split="_"))[1]))
    
    for(i in 1:length(files)){
      dataset=datasets[i]
      
      print(dataset)

      if(dataset!="VISTA"){
        
        path=paste(pathEnhancers, sp, "/coord_enh/", files[i], sep="")
        
        if(file.exists(path)){
          header=readLines(path, n=1)
          firstcol=unlist(strsplit(header, split="\t"))[1]
          
          if(firstcol=="chr"){
            this.coords=read.table(path, h=T, stringsAsFactors=F, sep="\t")
          } else{
            this.coords=read.table(path, h=F, stringsAsFactors=F, sep="\t")
            colnames(this.coords)[1]="chr"
            colnames(this.coords)[2]="start"
            colnames(this.coords)[3]="end"
          }
          
          this.coords$chr=unlist(lapply(this.coords$chr, function(x) substr(x, 4, nchar(x))))
          this.coords=this.coords[, c("chr", "start", "end")]
          
          enhancer.coords[[sp]][[dataset]]=this.coords
        }
      }
    }
  }
}

#####################################################################

save(enhancer.coords, file="RData/data.enhancers.RData")

#####################################################################

