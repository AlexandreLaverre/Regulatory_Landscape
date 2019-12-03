####################################################################################

exp.human=read.table("../../results/expression_estimation/Human/AllSamples_KallistoNormalizedTPM_FilteredTranscripts.txt", h=T, stringsAsFactors=F)
exp.mouse=read.table("../../results/expression_estimation/Mouse/AllSamples_KallistoNormalizedTPM_FilteredTranscripts.txt", h=T, stringsAsFactors=F)

rownames(exp.human)=exp.human$GeneID
rownames(exp.mouse)=exp.mouse$GeneID

exp.human=exp.human[,-1]
exp.mouse=exp.mouse[,-1]

colnames(exp.human)=paste("Human", colnames(exp.human), sep="_")
colnames(exp.mouse)=paste("Mouse", colnames(exp.mouse), sep="_")

####################################################################################

ortho=read.table("../../data/ensembl_ortho/Ortho_Human_Mouse_Ensembl94.txt", h=T, stringsAsFactors=F, sep="\t", quote="\"")
colnames(ortho)=c("Human", "Mouse", "HomologyType", "dN", "dS")
ortho=ortho[which(ortho$HomologyType=="ortholog_one2one" & !is.na(ortho$dN) & !is.na(ortho$dS)),]

ortho=ortho[which(ortho$Human%in%rownames(exp.human) & ortho$Mouse%in%rownames(exp.mouse)),]

####################################################################################

exp=cbind(exp.human[ortho$Human, ], exp.mouse[ortho$Mouse,], stringsAsFactors=F)

species=unlist(lapply(colnames(exp), function(x) unlist(strsplit(x, split="_"))[1]))
celltype=rep(NA, dim(exp)[2])
celltype[grep("B", colnames(exp))]="Bcell"
celltype[grep("ESC", colnames(exp))]="ESC"

pch.sp=c(21, 22)
names(pch.sp)=c("Human", "Mouse")

col.celltype=c("indianred", "steelblue")
names(col.celltype)=c("Bcell", "ESC")

####################################################################################

library(ade4)
library(ape)

####################################################################################

pca=dudi.pca(t(log2(exp+1)), center=T, scale=F, scannf=F, nf=5)

####################################################################################

plot(pca$li[,1], pca$li[,2], pch=pch.sp[species], col=col.celltype[celltype], bg=col.celltype[celltype])

####################################################################################

cormat=cor(exp, method="spearman")

tree=nj(1-cormat)
##tree=root(tree, outgroup=grep("ESC", tree$tip.label))

plot(tree)

####################################################################################

## average expression across species and replicates

fac.spcell=as.factor(paste(species, celltype, sep="_"))

exp.avg=as.data.frame(t(apply(exp, 1, function(x) tapply(as.numeric(x), fac.spcell, mean))))
exp.median=as.data.frame(t(apply(exp, 1, function(x) tapply(as.numeric(x), fac.spcell, median))))

####################################################################################
