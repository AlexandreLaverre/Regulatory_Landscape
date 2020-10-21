####################################################################################

options(stringsAsFactors=F)

####################################################################################

## define paths

user=as.character(Sys.getenv()["USER"])

if(user=="laverre"){
  pathFinalData="/home/laverre/Manuscript/"
}

if(user=="necsulea"){
  pathFinalData="/beegfs/data/necsulea/RegulatoryLandscapesManuscript/"
  pathScripts="/beegfs/data/necsulea/RegulatoryLandscapes/scripts/"
}

if(user=="ubuntu"){
  pathFinalData="/mnt/RegulatoryLandscapesManuscript/"
}

####################################################################################

pathExpressionData=paste(pathFinalData, "SupplementaryDataset6/", sep="")
pathOrtho=paste(pathFinalData, "SupplementaryDataset7/", sep="")

pathScriptsExpression=paste(pathScripts, "expression_estimation/", sep="")

####################################################################################

source(paste(pathScriptsExpression, "normalization.R", sep=""))

####################################################################################

exp.human=read.table(paste(pathExpressionData, "human/TPMValues_CommonCellTypes.txt", sep=""), h=T, stringsAsFactors=F)
exp.mouse=read.table(paste(pathExpressionData, "mouse/TPMValues_CommonCellTypes.txt", sep=""), h=T, stringsAsFactors=F)

rownames(exp.human)=exp.human$GeneID
rownames(exp.mouse)=exp.mouse$GeneID

exp.human=exp.human[,-which(colnames(exp.human)=="GeneID")]
exp.mouse=exp.mouse[,-which(colnames(exp.mouse)=="GeneID")]

colnames(exp.human)=paste("Human", colnames(exp.human), sep="_")
colnames(exp.mouse)=paste("Mouse", colnames(exp.mouse), sep="_")

####################################################################################

ortho=read.table(paste(pathOrtho, "human/gene_orthology/human2mouse_orthologue_dNdS.txt",sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
colnames(ortho)=c("Human", "Mouse", "HomologyType", "dN", "dS")
ortho=ortho[which(ortho$HomologyType=="ortholog_one2one" & !is.na(ortho$dN) & !is.na(ortho$dS)),]

ortho=ortho[which(ortho$Human%in%rownames(exp.human) & ortho$Mouse%in%rownames(exp.mouse)),]

####################################################################################

exp.ortho=cbind(exp.human[ortho$Human,], exp.mouse[ortho$Mouse,])
rownames(exp.ortho)=paste(ortho$Human, ortho$Mouse, sep="_")

norm.data=normalization(exp.ortho)
exp.ortho.norm=norm.data[["expdata.norm"]]
rownames(exp.ortho.norm)=rownames(exp.ortho.norm)

####################################################################################

samples=colnames(exp.ortho)

species=unlist(lapply(samples, function(x) unlist(strsplit(x, split="_"))[1]))

celltype=rep(NA, length(samples))
celltype[grep("adipo", samples)]="adipocytes"
celltype[grep("B", samples)]="Bcells"
celltype[grep("ESC", samples)]="ESC"

spcelltype=as.factor(paste(species, celltype, sep="_"))

####################################################################################

meanexp=t(apply(exp.ortho.norm, 1, function(x) tapply(as.numeric(x), spcelltype, mean)))
medianexp=t(apply(exp.ortho.norm, 1, function(x) tapply(as.numeric(x), spcelltype, median)))

####################################################################################

expdiv.median=list()
expdiv.mean=list()

for(cell in unique(celltype)){
  expdiv.median[[cell]]=abs(medianexp[,paste("Human",cell,sep="_")]-medianexp[,paste("Mouse",cell,sep="_")])/apply(medianexp[,paste(c("Human", "Mouse"), cell, sep="_")], 1, max)
  expdiv.mean[[cell]]=abs(meanexp[,paste("Human",cell,sep="_")]-meanexp[,paste("Mouse",cell,sep="_")])/apply(meanexp[,paste(c("Human", "Mouse"), cell, sep="_")], 1, max)
}

expdiv.median=as.data.frame(expdiv.median)
expdiv.mean=as.data.frame(expdiv.mean)

####################################################################################

write.table(expdiv.median, file=paste(pathExpressionData,"expression_divergence/ExpressionDivergence_CellTypes_MedianTPM.txt", sep=""), row.names=T, col.names=T, sep="\t", quote=F)
write.table(expdiv.mean, file=paste(pathExpressionData, "/expression_divergence/ExpressionDivergence_CellTypes_MeanTPM.txt", sep=""), row.names=T, col.names=T, sep="\t", quote=F)

####################################################################################
