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

ensrelease=94

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

####################################################################################

samples=colnames(exp.ortho)
species=unlist(lapply(samples, function(x) unlist(strsplit(x, split="_"))[1]))
celltype=rep(NA, length(samples))
celltype[grep("adipo", samples)]="adipocytes"
celltype[grep("B", samples)]="Bcells"
celltype[grep("ESC", samples)]="ESC"

pch.sp=c(21, 24)
names(pch.sp)=c("Human", "Mouse")

col.celltype=c("darkorange", "steelblue", "darkgreen")
names(col.celltype)=c("adipocytes", "Bcells", "ESC")

####################################################################################
  
library(ade4)

####################################################################################

pca=dudi.pca(log2(exp.ortho+1), center=T, scale=T, scannf=F, nf=5)

####################################################################################

explained=round(100*pca$eig/sum(pca$eig))

pdf(file="figures/PCA_CellTypes_HumanMouse.pdf", width=7, height=7)

par(mar=c(4.1, 4.1, 2.1, 1.1))
plot(pca$co[,1], pca$co[,2],  col="black", bg=col.celltype[celltype], pch=pch.sp[species], xlab=paste("coordinates on PC1 (", explained[1], "% explained variance)", sep=""),  ylab=paste("coordinates on PC2 (", explained[2], "% explained variance)", sep=""), cex=1.25)

legend("topright", legend=c("adipocytes", "B cells", "ESC"), fill=col.celltype, inset=0.01, bty="n")
legend("top", legend=c("human", "mouse"), pch=pch.sp, inset=c(0.2, 0.01), bty="n")

dev.off()

####################################################################################
