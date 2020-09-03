###########################################################################

source("parameters.R")

pathExpression=paste(pathFinalData, "SupplementaryDataset6/", sep="")

###########################################################################

exp.common.celltypes=list()

for(sp in c("human", "mouse")){
  this.exp=read.table(paste(pathExpression, sp, "/TPMValues_CommonCellTypes.txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
  rownames(this.exp)=this.exp$GeneID
   
  this.exp=this.exp[,which(colnames(this.exp)!="GeneID")]
  this.exp=as.matrix(this.exp)

  
  celltype=rep(NA, dim(this.exp)[2])
  celltype[grep("adipo", colnames(this.exp))]="adipo"
  celltype[grep("ESC", colnames(this.exp))]="ESC"
  celltype[grep("B", colnames(this.exp))]="B"
  celltype=as.factor(celltype)
  
  avg.exp=t(apply(this.exp, 1, function(x) tapply(x, celltype, mean)))
  
  exp.common.celltypes[[sp]]=list("allsamples"=this.exp, "avg"=avg.exp)
}

###########################################################################

expstats.cm2019=list()

for(sp in c("human", "mouse")){
  this.exp=read.table(paste(pathExpression, sp, "/ExpressionStatistics_CardosoMoreira2019.txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
  rownames(this.exp)=this.exp$GeneID

  expstats.cm2019[[sp]]=this.exp
}

###########################################################################

avgexp.cm2019=list()

for(sp in c("human", "mouse")){
  this.exp=read.table(paste(pathExpression, sp, "/AverageRPKM_CardosoMoreira2019.txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
 
  avgexp.cm2019[[sp]]=this.exp
}

###########################################################################

save(list=c("exp.common.celltypes", "expstats.cm2019", "avgexp.cm2019"), file=paste(pathFigures, "RData/data.gene.expression.RData", sep=""))

###########################################################################
