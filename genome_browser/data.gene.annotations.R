#####################################################################

pathAnouk="/beegfs/data/necsulea/RegulatoryLandscapes_AN/"
pathAlex="/beegfs/data/alaverre/Regulatory_landscape/"

pathAnnotations=paste(pathAnouk, "data/ensembl_annotations/", sep="")

ensrelease=94

options(stringsAsFactors=F)

#####################################################################

gene.coords=list()
exon.coords=list()

for(sp in c("Human", "Mouse")){
  print(sp)
  
  ## gene coordinates
  print("gene coordinates")
  
  annot=read.table(paste(pathAnnotations, sp, "/GeneCoordinates_Ensembl94.gtf", sep=""), h=F, sep="\t")
  colnames(annot)=c("chr", "source", "type", "start", "end", "frame", "strand", "score", "info")

  info=lapply(annot$info, function(x) unlist(strsplit(x, split=";")))
  geneid=unlist(lapply(info, function(x) {y=grep("gene_id", x, value=T); if(length(y)==1){return(unlist(strsplit(y[1], split=" "))[[2]])} else{return(NA)}}))
  biotype=unlist(lapply(info, function(x) {y=grep("gene_biotype", x, value=T); if(length(y)==1){return(unlist(strsplit(y[1], split=" "))[[3]])} else{return(NA)}}))
  name=unlist(lapply(info, function(x) {y=grep("gene_name", x, value=T); if(length(y)==1){return(unlist(strsplit(y[1], split=" "))[[3]])} else{return(NA)}}))

  annot$id=geneid
  annot$biotype=biotype
  annot$name=name
  
  gene.coords[[tolower(sp)]]=annot

  ## exon coordinates
  print("exon coordinates")

  annot.exons=read.table(paste(pathAnnotations, sp, "/ExonCoordinates_Ensembl94.gtf", sep=""), h=F, sep="\t")
  colnames(annot.exons)=c("chr", "source", "type", "start", "end", "frame", "strand", "score", "info")

  info=lapply(annot.exons$info, function(x) unlist(strsplit(x, split=";")))
  geneid=unlist(lapply(info, function(x) {y=grep("gene_id", x, value=T); if(length(y)==1){return(unlist(strsplit(y[1], split=" "))[[2]])} else{return(NA)}}))

  annot.exons$geneid=geneid

  annot.exons=annot.exons[,c("geneid", "chr", "start", "end", "strand")]
  
  exon.coords[[tolower(sp)]]=annot.exons
}

#####################################################################

save(list=c("gene.coords", "exon.coords"), file="RData/data.gene.annotations.RData")

#####################################################################

