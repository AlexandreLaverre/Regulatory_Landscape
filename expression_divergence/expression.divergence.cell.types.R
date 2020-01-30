####################################################################################

options(stringsAsFactors=F)

####################################################################################

pwd=getwd()
dirs=unlist(strsplit(pwd, split="/"))
path=paste(dirs[1:(length(dirs)-2)], collapse="/")
path=paste(path, "/", sep="")

####################################################################################

pathExpression=paste(path, "results/expression_estimation/", sep="")
pathOrtho=paste(path, "data/ensembl_ortho/", sep="")
pathScriptsExpression=paste(path, "scripts/expression_estimation/", sep="")

ensrelease=94

####################################################################################

source(paste(pathScriptsExpression, "normalization.R", sep=""))

####################################################################################

exp.human=read.table(paste(pathExpression, "Human/AllSamples_KallistoNormalizedTPM_FilteredTranscripts.txt", sep=""), h=T, stringsAsFactors=F)
exp.mouse=read.table(paste(pathExpression, "Mouse/AllSamples_KallistoNormalizedTPM_FilteredTranscripts.txt", sep=""), h=T, stringsAsFactors=F)

rownames(exp.human)=exp.human$GeneID
rownames(exp.mouse)=exp.mouse$GeneID

exp.human=exp.human[,-1]
exp.mouse=exp.mouse[,-1]

colnames(exp.human)=paste("Human", colnames(exp.human), sep="_")
colnames(exp.mouse)=paste("Mouse", colnames(exp.mouse), sep="_")

####################################################################################

ortho=read.table(paste(pathOrtho, "Ortho_Human_Mouse_Ensembl",ensrelease,".txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
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

write.table(exp.ortho, file=paste(pathExpression, "AllSpecies_AllSamples_OrthoGenes_KallistoTPM_FilteredTranscripts.txt", sep=""), row.names=T, col.names=T, sep="\t")

write.table(exp.ortho.norm, file=paste(pathExpression, "AllSpecies_AllSamples_OrthoGenes_KallistoTPM_FilteredTranscripts_NormalizedHK.txt", sep=""), row.names=T, col.names=T, sep="\t")
            

####################################################################################

samples=colnames(exp.ortho)

species=unlist(lapply(samples, function(x) unlist(strsplit(x, split="_"))[1]))

celltype=rep(NA, length(samples))
celltype[grep("adipo", samples)]="adipocytes"
celltype[grep("B", samples)]="Bcells"
celltype[grep("ESC", samples)]="ESC"

####################################################################################



####################################################################################
  

