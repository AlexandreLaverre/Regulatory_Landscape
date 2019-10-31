#####################################################################

pathAnouk="/beegfs/data/necsulea/RegulatoryLandscapes/"
pathAlex="/beegfs/data/alaverre/Regulatory_landscape/"

pathAnnotations=paste(pathAnouk, "data/ensembl_annotations/", sep="")

ensrelease=94

options(stringsAsFactors=F)

#####################################################################

gene.coords=list()

for(sp in c("Human", "Mouse")){
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
}

#####################################################################

save(gene.coords, file="RData/data.annotations.RData")

#####################################################################

