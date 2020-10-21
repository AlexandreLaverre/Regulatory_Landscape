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

######################################################################################

## we take into account only somatic organs

samples=kronecker(c("Brain", "Cerebellum", "Kidney", "Liver", "Heart"), unique(stages$ID), paste, sep="_")

samples=samples[which(samples%in%colnames(avgexp.human) & samples%in%colnames(avgexp.mouse))]

avgexp.human=avgexp.human[,samples]
avgexp.mouse=avgexp.mouse[,samples]

######################################################################################

## only ortho genes

avgexp.human=avgexp.human[ortho$Human,]
avgexp.mouse=avgexp.mouse[ortho$Mouse,]

######################################################################################

## select genes that are expressed above an RPKM threshold in at least some samples

minRPKM=2
minsamples=3

nbsamples.human=apply(avgexp.human, 1, function(x) length(which(x>=minRPKM)))
nbsamples.mouse=apply(avgexp.mouse, 1, function(x) length(which(x>=minRPKM)))

okgenes=which(nbsamples.human>=minsamples & nbsamples.mouse>=minsamples)

avgexp.human=avgexp.human[okgenes,]
avgexp.mouse=avgexp.mouse[okgenes,]

######################################################################################

relexp.human=avgexp.human/apply(avgexp.human,1,sum)
relexp.mouse=avgexp.mouse/apply(avgexp.mouse,1,sum)

relexp.human=as.matrix(relexp.human)
relexp.mouse=as.matrix(relexp.mouse)

######################################################################################

distance=unlist(lapply(1:dim(relexp.human)[1], function(x) {y=as.numeric(relexp.mouse[x,]); z=as.numeric(relexp.human[x,]); return(sqrt(sum((y-z)^2)))}))
names(distance)=rownames(relexp.human)

######################################################################################

correlation.spearman=unlist(lapply(1:dim(relexp.human)[1], function(x) cor(relexp.human[x,], relexp.mouse[x,], method="spearman")))
names(correlation.spearman)=rownames(relexp.human)

correlation.pearson=unlist(lapply(1:dim(relexp.human)[1], function(x) cor(relexp.human[x,], relexp.mouse[x,], method="pearson")))
names(correlation.pearson)=rownames(relexp.human)

######################################################################################

results=data.frame("IDMouse"=rownames(relexp.mouse), "IDHuman"=rownames(relexp.human), "EuclideanDistance"=distance, "CorrelationSpearman"=correlation.spearman, "CorrelationPearson"=correlation.pearson, stringsAsFactors=F)

######################################################################################

## correct for average expression 

for(sample in colnames(avgexp.human)){
  results[,paste("Human_",sample,sep="")]=avgexp.human[,sample]
}

for(sample in colnames(avgexp.mouse)){
  results[,paste("Mouse_",sample,sep="")]=avgexp.mouse[,sample]
}

results$Human_MeanRPKM=apply(results[,grep("^Human_", colnames(results))],1, mean)
results$Mouse_MeanRPKM=apply(results[,grep("^Mouse_", colnames(results))],1, mean)

results$MeanRPKM=(results$Human_MeanRPKM+results$Mouse_MeanRPKM)/2

lm1=lm(results$EuclideanDistance~log2(results$MeanRPKM+1))
results$ResidualEuclideanDistance=lm1$residuals

lm2=lm(results$CorrelationSpearman~log2(results$MeanRPKM+1))
results$ResidualSpearman=lm2$residuals

lm3=lm(results$CorrelationPearson~log2(results$MeanRPKM+1))
results$ResidualPearson=lm3$residuals

######################################################################################

write.table(results, file=paste(paste(pathSupplementaryDataset6, "expression_divergence/ExpressionDivergence_CardosoMoreira2019.txt", row.names=F, col.names=T, sep="\t", quote=F)

######################################################################################
######################################################################################
