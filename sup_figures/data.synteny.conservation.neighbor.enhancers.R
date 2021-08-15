##########################################################################

options(stringsAsFactors = FALSE)

#########################################################################

## if it's the first time we run this figure, we load and prepare data

objects=ls()

if(!"pathScripts"%in%objects){
  

  source("../main_figures/parameters.R")
  
  minAlnLength=10

  set.seed(19)
}

##################################################################

load(paste(pathFigures, "RData/data.gene.annotations.RData",sep=""))
load(paste(pathFigures, "RData/data.neighbor.enhancers.RData",sep=""))
load(paste(pathFigures, "RData/data.ortho.genes.RData", sep=""))

##################################################################

enh="ENCODE"

## target chromosomes for enhancers

enhancer.chr=list()

for(ref in c("human", "mouse")){
  tg=setdiff(c("human", "mouse"), ref)
  
  this.aln=read.table(paste(pathFinalData, "SupplementaryDataset7/", ref,"/sequence_conservation/enhancers/ENCODE/AlignmentStatistics_Excluding_Exons_",ref,"2",tg, ".txt", sep=""), h=T, sep="\t", stringsAsFactors=F)
  
  rownames(this.aln)=this.aln[,paste("ID.", ref, sep="")]
  
  ## filter on sequence alignment
  this.aln$pcungapped=this.aln$FilteredUngappedLength/this.aln$FilteredAlignmentLength
  this.aln$pcungapped[which(this.aln$FilteredAlignmentLength < minAlnLength)]=NA

  ## 10% quantile
  align.threshold=quantile(this.aln$pcungapped, p=0.1, na.rm=T)
  this.aln=this.aln[which(this.aln$pcungapped>=align.threshold),]

  ## chromosome in target species
  this.aln[,paste("chr",ref, sep=".")]=unlist(lapply(this.aln[,paste("ID.", ref, sep="")], function(x) unlist(strsplit(x, split=":"))[1]))
  this.aln[,paste("chr",tg, sep=".")]=unlist(lapply(this.aln[,paste("ID.", tg, sep="")], function(x) unlist(strsplit(x, split=":"))[1]))

  this.chr=this.aln[,paste("chr", tg, sep=".")]
  names(this.chr)=rownames(this.aln)

  ## remove the "^chr

  g=grep("^chr", this.chr)
  if(length(g)>0){
    this.chr=unlist(lapply(this.chr, function(x) unlist(strsplit(x, split="chr"))[2]))
  }
  
  enhancer.chr[[paste(ref, "2", tg, sep="")]]=this.chr
}

##################################################################

## target chromosomes for genes

gene.chr=list()

for(ref in c("human", "mouse")){
  tg=setdiff(c("human", "mouse"), ref)

  this.annot=gene.annot[[tg]]

  this.chr=this.annot[ortho[,tg],"Chr"]
  names(this.chr)=ortho[,ref]

  gene.chr[[paste(ref, "2", tg, sep="")]]=this.chr
}

##################################################################

## check synteny conservation

synteny.cons=list()

for(ref in c("human", "mouse")){
  tg=setdiff(c("human", "mouse"), ref)

  contacts=gene.enhancer.contacts[[ref]][["ENCODE"]][["real"]]

  this.gene.chr=gene.chr[[paste(ref, "2", tg, sep="")]]
  this.enhancer.chr=enhancer.chr[[paste(ref, "2", tg, sep="")]]
  
  contacts[,"TargetChrGene"]=this.gene.chr[contacts$GeneID]
  contacts[,"TargetChrEnhancer"]=this.enhancer.chr[contacts$EnhancerID]

  contacts[,"ConservedSynteny"]=rep(NA, dim(contacts)[1])
  contacts[which(contacts[,"TargetChrGene"]==contacts[,"TargetChrEnhancer"]),"ConservedSynteny"]="yes"
  contacts[which(contacts[,"TargetChrGene"]!=contacts[,"TargetChrEnhancer"]),"ConservedSynteny"]="no"

  contacts=contacts[which(!is.na(contacts$ConservedSynteny)),]
  
  synteny.cons[[paste(ref, "2", tg, sep="")]]=contacts
}

##################################################################

save(synteny.cons, file=paste(pathFigures, "RData/data.synteny.conservation.neighbor.enhancers.RData",sep=""))

##################################################################
