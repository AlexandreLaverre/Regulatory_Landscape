###########################################################################

source("parameters.R")

pathEvolution=paste(pathFinalData, "SupplementaryDataset7/", sep="")

###########################################################################
## only for human and mouse

all.ortho=read.table(paste(pathEvolution, "human/gene_orthology/human2mouse_orthologue.txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")

###########################################################################

ortho=all.ortho[which(all.ortho$Mouse.homology.type=="ortholog_one2one"), c("Gene.stable.ID", "Mouse.gene.stable.ID")]
colnames(ortho)=c("human", "mouse")

###########################################################################

save(list=c("ortho"), file=paste(pathFigures, "RData/data.ortho.genes.RData", sep=""))

###########################################################################
