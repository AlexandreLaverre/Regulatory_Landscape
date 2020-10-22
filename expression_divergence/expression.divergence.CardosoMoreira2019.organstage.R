######################################################################################

## define paths

user=as.character(Sys.getenv()["USER"])

if(user=="laverre"){
  pathFinalData="/home/laverre/Manuscript/"
}

if(user=="necsulea"){
  pathFinalData="/beegfs/data/necsulea/RegulatoryLandscapesManuscript/"
}


if(user=="ubuntu"){
  pathFinalData="/mnt/RegulatoryLandscapesManuscript/"
}

######################################################################################

pathExpressionData=paste(pathFinalData, "SupplementaryDataset6/", sep="")
pathOrtho=paste(pathFinalData, "SupplementaryDataset7/", sep="")
 
minRPKM=2
minsamples=3

######################################################################################

## correspondence between developmental stages
stages=read.table(paste(pathExpressionData, "expression_divergence/DevelopmentalStageCorrespondence_CardosoMoreira2019.txt", sep=""), h=T, stringsAsFactors=F)

######################################################################################

exp.human=read.table(paste(pathExpressionData, "human/AverageRPKM_CardosoMoreira2019.txt", sep=""), h=T, stringsAsFactors=F)
exp.mouse=read.table(paste(pathExpressionData, "mouse/AverageRPKM_CardosoMoreira2019.txt", sep=""), h=T, stringsAsFactors=F)
  
######################################################################################

ortho=read.table(paste(pathOrtho, "human/gene_orthology/human2mouse_orthologue_dNdS.txt",sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
colnames(ortho)=c("Human", "Mouse", "HomologyType", "dN", "dS")
ortho=ortho[which(ortho$HomologyType=="ortholog_one2one" & !is.na(ortho$dN) & !is.na(ortho$dS)),]

ortho=ortho[which(ortho$Human%in%rownames(exp.human) & ortho$Mouse%in%rownames(exp.mouse)),]

######################################################################################

stage.human=unlist(lapply(colnames(exp.human), function(x) unlist(strsplit(x, split="\\."))[2]))
stage.mouse=unlist(lapply(colnames(exp.mouse), function(x) {y=unlist(strsplit(x, split="\\.")); return(paste(y[-1], collapse="."))}))

  ## select only comparable stages

exp.human=exp.human[,which(stage.human%in%stages$Human)]
exp.mouse=exp.mouse[,which(stage.mouse%in%stages$Mouse)]

stage.human=stage.human[which(stage.human%in%stages$Human)]
stage.mouse=stage.mouse[which(stage.mouse%in%stages$Mouse)]

## change stage name to stage ID

orthostage.human=unlist(lapply(stage.human, function(x) paste(unique(stages$ID[which(stages$Human==x)]), collapse=",")))
orthostage.mouse=unlist(lapply(stage.mouse, function(x) paste(unique(stages$ID[which(stages$Mouse==x)]), collapse=",")))

######################################################################################

organ.human=unlist(lapply(colnames(exp.human), function(x) unlist(strsplit(x, split="\\."))[1]))
organ.mouse=unlist(lapply(colnames(exp.mouse), function(x) unlist(strsplit(x, split="\\."))[1]))

organstage.human=paste(organ.human, orthostage.human, sep="_")
organstage.mouse=paste(organ.mouse, orthostage.mouse, sep="_")

## select common organ/stage combinations

common=intersect(organstage.human, organstage.mouse)

exp.human=exp.human[,which(organstage.human%in%common)]
exp.mouse=exp.mouse[,which(organstage.mouse%in%common)]

organstage.human=organstage.human[which(organstage.human%in%common)]
organstage.mouse=organstage.mouse[which(organstage.mouse%in%common)]

######################################################################################

exp.human=as.matrix(exp.human)
exp.mouse=as.matrix(exp.mouse)

organstage.human=as.factor(organstage.human)
organstage.mouse=as.factor(organstage.mouse)

avgexp.human=t(apply(exp.human, 1, function(x) tapply(x, organstage.human, mean)))
avgexp.mouse=t(apply(exp.mouse, 1, function(x) tapply(x, organstage.mouse, mean)))

medexp.human=t(apply(exp.human, 1, function(x) tapply(x, organstage.human, median)))
medexp.mouse=t(apply(exp.mouse, 1, function(x) tapply(x, organstage.mouse, median)))

######################################################################################

samples=kronecker(c("Brain", "Cerebellum", "Kidney", "Liver", "Heart", "Ovary", "Testis"), unique(stages$ID), paste, sep="_")

samples=samples[which(samples%in%colnames(avgexp.human) & samples%in%colnames(avgexp.mouse))]

avgexp.human=avgexp.human[,samples]
avgexp.mouse=avgexp.mouse[,samples]

medexp.human=medexp.human[,samples]
medexp.mouse=medexp.mouse[,samples]

######################################################################################

## only ortho genes

avgexp.human=avgexp.human[ortho$Human,]
avgexp.mouse=avgexp.mouse[ortho$Mouse,]

medexp.human=medexp.human[ortho$Human,]
medexp.mouse=medexp.mouse[ortho$Mouse,]

######################################################################################

expdiv.median=list()
expdiv.mean=list()

for(sample in unique(samples)){
  expdiv.median[[sample]]=abs(medexp.human[,sample]-medexp.mouse[,sample])/apply(cbind(medexp.human[,sample], medexp.mouse[,sample]), 1, max)
  expdiv.mean[[sample]]=abs(avgexp.human[,sample]-avgexp.mouse[,sample])/apply(cbind(avgexp.human[,sample], avgexp.mouse[,sample]), 1, max)
}

expdiv.median=as.data.frame(expdiv.median)
expdiv.mean=as.data.frame(expdiv.mean)

######################################################################################

write.table(expdiv.median, file=paste(pathExpressionData,"expression_divergence/ExpressionDivergence_CardosoMoreira2019_OrganStage_MedianTPM.txt", sep=""), row.names=T, col.names=T, sep="\t", quote=F)

write.table(expdiv.mean, file=paste(pathExpressionData, "/expression_divergence/ExpressionDivergence_CardosoMoreira2019_OrganStage__MeanTPM.txt", sep=""), row.names=T, col.names=T, sep="\t", quote=F)

######################################################################################
######################################################################################
