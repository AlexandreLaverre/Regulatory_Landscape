######################################################################################
#setwd("/home/laverre/Anouk/data/enhancers/")

stages=read.table("../../data/gene_expression_CardosoMoreira2019/DevelopmentalStages.txt", h=T, stringsAsFactors=F)
exp.human=read.table("../../data/gene_expression_CardosoMoreira2019/Human.RPKM.txt", h=T, stringsAsFactors=F)
exp.mouse=read.table("../../data/gene_expression_CardosoMoreira2019/Mouse.RPKM.txt", h=T, stringsAsFactors=F)

######################################################################################

ortho=read.table("../../data/ensembl_ortho/Ortho_Human_Mouse_Ensembl94.txt", h=T, stringsAsFactors=F, sep="\t", quote="\"")
colnames(ortho)=c("Human", "Mouse", "HomologyType", "dN", "dS")
ortho=ortho[which(ortho$HomologyType=="ortholog_one2one" & !is.na(ortho$dN) & !is.na(ortho$dS)),]

######################################################################################

stage.human=unlist(lapply(colnames(exp.human), function(x) unlist(strsplit(x, split="\\."))[2]))
stage.mouse=unlist(lapply(colnames(exp.mouse), function(x) {y=unlist(strsplit(x, split="\\.")); n=length(y); return(paste(y[2:(n-1)], collapse="."))}))

organ.human=unlist(lapply(colnames(exp.human), function(x) unlist(strsplit(x, split="\\."))[1]))
organ.mouse=unlist(lapply(colnames(exp.mouse), function(x) unlist(strsplit(x, split="\\."))[1]))

organstage.human=as.factor(paste(organ.human, stage.human, sep="_"))
organstage.mouse=as.factor(paste(organ.mouse, stage.mouse, sep="_"))

exp.human=as.matrix(exp.human)
exp.mouse=as.matrix(exp.mouse)

######################################################################################

avgexp.human=t(apply(exp.human, 1, function(x) tapply(x, organstage.human, mean)))
avgexp.mouse=t(apply(exp.mouse, 1, function(x) tapply(x, organstage.mouse, mean)))

######################################################################################

samples.human=kronecker(c("Brain", "Cerebellum", "Kidney", "Liver", "Heart"), stages$Human, paste, sep="_")
samples.mouse=kronecker(c("Brain", "Cerebellum", "Kidney", "Liver", "Heart"), stages$Mouse, paste, sep="_")

ok=which(samples.human%in%colnames(avgexp.human) & samples.mouse%in%colnames(avgexp.mouse))
samples.human=samples.human[ok]
samples.mouse=samples.mouse[ok]

avgexp.human=avgexp.human[,samples.human]
avgexp.mouse=avgexp.mouse[,samples.mouse]

######################################################################################

ortho=ortho[which(ortho$Human%in%rownames(avgexp.human) & ortho$Mouse%in%rownames(avgexp.mouse)),]

avgexp.human=avgexp.human[ortho$Human,]
avgexp.mouse=avgexp.mouse[ortho$Mouse,]

######################################################################################

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

######################################################################################

distance=unlist(lapply(1:dim(relexp.human)[1], function(x) {y=as.numeric(relexp.mouse[x,]); z=as.numeric(relexp.human[x,]); return(sqrt(sum((y-z)^2)))}))
names(distance)=rownames(relexp.human)

######################################################################################

correlation.sper=unlist(lapply(1:dim(relexp.human)[1],
                          function(x) {y=as.numeric(relexp.mouse[x,]); z=as.numeric(relexp.human[x,]); return(abs(cor(y,z, method="spearman")))}))
names(correlation.sper)=rownames(relexp.human)

correlation.pear=unlist(lapply(1:dim(relexp.human)[1],
                               function(x) {y=as.numeric(relexp.mouse[x,]); z=as.numeric(relexp.human[x,]); return(abs(cor(y,z, method="pearson")))}))
names(correlation.pear)=rownames(relexp.human)

######################################################################################

results=data.frame("IDMouse"=rownames(relexp.mouse), "IDHuman"=rownames(relexp.human), "ExpressionDivergence"=distance,
                   "CorrelationSpearman"=correlation.sper, "CorrelationPearson"=correlation.pear, stringsAsFactors=F)


for(sample in colnames(avgexp.human)){
  results[,paste("Human_",sample,sep="")]=avgexp.human[,sample]
}

for(sample in colnames(avgexp.mouse)){
  results[,paste("Mouse_",sample,sep="")]=avgexp.mouse[,sample]
}

results$Human_MeanRPKM=apply(results[,grep("^Human_", colnames(results))],1, mean)
results$Mouse_MeanRPKM=apply(results[,grep("^Mouse_", colnames(results))],1, mean)

results$MeanRPKM=(results$Human_MeanRPKM+results$Mouse_MeanRPKM)/2

lm1=lm(results$ExpressionDivergence~log2(results$MeanRPKM+1))
results$ResidualExpressionDivergence=lm1$residuals

lm1=lm(results$CorrelationSpearman~log2(results$MeanRPKM+1))
results$ResidualSpearman=lm1$residuals

lm1=lm(results$CorrelationPearson~log2(results$MeanRPKM+1))
results$ResidualPearson=lm1$residuals

######################################################################################

write.table(results, file="/home/laverre/Manuscript/SupplementaryDataset6/expression_divergence/ExpressionDivergence_CardosoMoreira2019_correlations.txt", row.names=F, col.names=T, sep="\t", quote=F)

######################################################################################
######################################################################################
