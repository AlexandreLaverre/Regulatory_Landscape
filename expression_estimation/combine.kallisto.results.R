########################################################################

path="/mnt/RegulatoryLandscapes/"
pathExpression=paste(path,"results/expression_estimation/", sep="")
pathEnsembl=paste(path,"data/ensembl_annotations/", sep="")
pathStringTie=paste(path,"results/stringtie_assembly/",sep="")

set.seed(19)

source("normalization.R")

types=c("AllTranscripts", "FilteredTranscripts")
splist=c("Mouse", "Human")

########################################################################

for(sp in splist){
  print(sp)
  
  samples=system(paste("ls ", pathExpression, sp, " | grep -v txt", sep=""), intern=T)

  for(type in types){
    print(type)

    read.counts=list()
    tpm=list()
    
    for(sample in samples){

      print(sample)

      this.data=c()
      
      try(this.data<-read.table(paste(pathExpression, sp, "/", sample,"/kallisto_",type, "/abundance.tsv",sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\""))

      if(length(this.data)>0){
        this.data$GeneID=unlist(lapply(this.data$target_id, function(x) unlist(strsplit(x, split=":"))[1]))
      
        gene.counts=as.numeric(tapply(as.numeric(this.data$est_counts), as.factor(this.data$GeneID), sum, na.rm=T))
        gene.tpm=as.numeric(tapply(as.numeric(this.data$tpm), as.factor(this.data$GeneID), sum, na.rm=T))

        names(gene.counts)=levels(as.factor(this.data$GeneID))
        names(gene.tpm)=levels(as.factor(this.data$GeneID))

        if(max(gene.tpm)!=0){
          read.counts[[sample]]=gene.counts
          tpm[[sample]]=gene.tpm
        }
      }
    }

    samples=names(read.counts)
    
    ## reorder values

    gene.order=names(read.counts[[samples[1]]])
    
    for(sample in samples){
      print(paste("reordering", sample))
      read.counts[[sample]]=read.counts[[sample]][gene.order]
      tpm[[sample]]=tpm[[sample]][gene.order]
    }
    
    ## make data frames
    
    read.counts=as.data.frame(read.counts)
    rownames(read.counts)=gene.order
    
    tpm=as.data.frame(tpm)
    rownames(tpm)=gene.order

    norm.data=normalization(tpm)
    tpm.norm=norm.data[["expdata.norm"]]
    rownames(tpm.norm)=gene.order

    hk.genes=norm.data[["hk.genes"]]
    
    ## add gene id as a column
    
    tpm$GeneID=gene.order
    tpm.norm$GeneID=gene.order
    read.counts$GeneID=gene.order

    tpm=tpm[,c("GeneID", setdiff(colnames(tpm), "GeneID"))]
    tpm.norm=tpm.norm[,c("GeneID", setdiff(colnames(tpm), "GeneID"))]
    read.counts=read.counts[,c("GeneID", setdiff(colnames(tpm), "GeneID"))]

    ## write output 
    
    write.table(read.counts, file=paste(pathExpression, sp, "/AllSamples_KallistoEstimatedCounts_", type,".txt",sep=""), row.names=F, col.names=T, quote=F, sep="\t")
    
    write.table(tpm, file=paste(pathExpression, sp, "/AllSamples_KallistoTPM_", type,".txt",sep=""), row.names=F, col.names=T, quote=F, sep="\t")

    write.table(tpm.norm, file=paste(pathExpression, sp, "/AllSamples_KallistoNormalizedTPM_", type,".txt",sep=""), row.names=F, col.names=T, quote=F, sep="\t")

 
     
  }
}

########################################################################
