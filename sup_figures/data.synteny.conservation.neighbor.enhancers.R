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
enhancer.pos=list()

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
  this.aln[,paste("chr",tg, sep=".")]=unlist(lapply(this.aln[,paste("ID.", tg, sep="")], function(x) unlist(strsplit(x, split=":"))[1]))

  this.chr=this.aln[,paste("chr", tg, sep=".")]
  names(this.chr)=rownames(this.aln)

  ## remove the "^chr

  g=grep("^chr", this.chr)
  if(length(g)>0){
    this.chr=unlist(lapply(this.chr, function(x) unlist(strsplit(x, split="chr"))[2]))
  }
  
  ## position in target species
  this.aln[,paste("start",tg, sep=".")]=unlist(lapply(this.aln[,paste("ID.", tg, sep="")], function(x) as.numeric(unlist(strsplit(x, split=":"))[2])))
  this.aln[,paste("end",tg, sep=".")]=unlist(lapply(this.aln[,paste("ID.", tg, sep="")], function(x) as.numeric(unlist(strsplit(x, split=":"))[3])))

  this.pos=(this.aln[,paste("start",tg,sep=".")]+this.aln[,paste("start",tg,sep=".")])/2
  names(this.pos)=rownames(this.aln)
  
  enhancer.chr[[paste(ref, "2", tg, sep="")]]=this.chr
  enhancer.pos[[paste(ref, "2", tg, sep="")]]=this.pos
}

##################################################################

## target chromosomes for genes

gene.chr=list()
gene.pos=list()

for(ref in c("human", "mouse")){
  tg=setdiff(c("human", "mouse"), ref)

  this.annot=gene.annot[[tg]]

  this.chr=this.annot[ortho[,tg],"Chr"]
  names(this.chr)=ortho[,ref]

  this.strand=this.annot[ortho[,tg],"Strand"]
  this.start=this.annot[ortho[,tg],"Start"]
  this.end=this.annot[ortho[,tg],"End"]

  this.tss=rep(NA, length(this.strand))
  this.tss[which(this.strand==1)]=this.start[which(this.strand==1)]
  this.tss[which(this.strand==-1)]=this.end[which(this.strand==-1)]

  names(this.tss)=ortho[,ref]
  
  gene.chr[[paste(ref, "2", tg, sep="")]]=this.chr
  gene.pos[[paste(ref, "2", tg, sep="")]]=this.tss
}

##################################################################

## check synteny conservation

synteny.cons=list()

for(ref in c("human", "mouse")){
  tg=setdiff(c("human", "mouse"), ref)

  contacts=gene.enhancer.contacts[[ref]][["ENCODE"]][["real"]]

  this.gene.chr=gene.chr[[paste(ref, "2", tg, sep="")]]
  this.gene.pos=gene.pos[[paste(ref, "2", tg, sep="")]]
  
  this.enhancer.chr=enhancer.chr[[paste(ref, "2", tg, sep="")]]
  this.enhancer.pos=enhancer.pos[[paste(ref, "2", tg, sep="")]]

  
  contacts[,"TargetChrGene"]=this.gene.chr[contacts$GeneID]
  contacts[,"TargetChrEnhancer"]=this.enhancer.chr[contacts$EnhancerID]
  contacts[,"TargetDistance"]=abs(this.gene.pos[contacts$GeneID]-this.enhancer.pos[contacts$EnhancerID])

  contacts[,"ConservedSynteny"]=rep(NA, dim(contacts)[1])
  contacts[which(contacts[,"TargetChrGene"]==contacts[,"TargetChrEnhancer"] & contacts[,"TargetDistance"]<=maxDistanceSyntenyTarget),"ConservedSynteny"]="yes"
  contacts[which(contacts[,"TargetChrGene"]!=contacts[,"TargetChrEnhancer"] | contacts[,"TargetDistance"]>maxDistanceSyntenyTarget),"ConservedSynteny"]="no"

  contacts=contacts[which(!is.na(contacts$ConservedSynteny)),]
  
  synteny.cons[[paste(ref, "2", tg, sep="")]]=contacts
}

##################################################################

save(synteny.cons, file=paste(pathFigures, "RData/data.synteny.conservation.neighbor.enhancers.RData",sep=""))

##################################################################
