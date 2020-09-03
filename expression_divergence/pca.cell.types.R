####################################################################################

options(stringsAsFactors=F)

####################################################################################

pwd=getwd()
dirs=unlist(strsplit(pwd, split="/"))
path=paste(dirs[1:(length(dirs)-2)], collapse="/")
path=paste(path, "/", sep="")

####################################################################################

pathExpression=paste(path, "results/expression_estimation/", sep="")

ensrelease=94

####################################################################################

exp.ortho=read.table(paste(pathExpression, "AllSpecies_AllSamples_OrthoGenes_KallistoTPM_FilteredTranscripts_NormalizedHK.txt", sep=""), h=T, sep="\t", stringsAsFactors=F)

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

pdf(paste(path, "scripts/main_figures/old/PCA_CellTypes_HumanMouse.pdf", sep=""), width=7, height=7)

par(mar=c(4.1, 4.1, 2.1, 1.1))
plot(pca$co[,1], pca$co[,2],  col="black", bg=col.celltype[celltype], pch=pch.sp[species], xlab=paste("coordinates on PC1 (", explained[1], "% explained variance)", sep=""),  ylab=paste("coordinates on PC2 (", explained[2], "% explained variance)", sep=""), cex=1.25)

legend("topright", legend=c("adipocytes", "B cells", "ESC"), fill=col.celltype, inset=0.01)
legend("bottomright", legend=c("human", "mouse"), pch=pch.sp, inset=0.01)

dev.off()

####################################################################################
