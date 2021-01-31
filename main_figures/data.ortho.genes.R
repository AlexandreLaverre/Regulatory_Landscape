###########################################################################

source("parameters.R")

pathEvolution=paste(pathFinalData, "SupplementaryDataset7/", sep="")

###########################################################################

## only for human and mouse

all.ortho=read.table(paste(pathEvolution, "human/gene_orthology/human2mouse_orthologue_dNdS.txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
colnames(all.ortho)=c("IDHuman", "IDMouse", "HomologyType", "dN", "dS")

all.ortho$dN=as.numeric(all.ortho$dN)
all.ortho$dS=as.numeric(all.ortho$dS)
all.ortho$dNdS=all.ortho$dN/all.ortho$dS

###########################################################################

ortho=all.ortho[which(all.ortho$HomologyType=="ortholog_one2one" & all.ortho$dNdS<1), c("IDHuman", "IDMouse")]
colnames(ortho)=c("human", "mouse")

###########################################################################

save(ortho, file=paste(pathFigures, "RData/data.ortho.genes.RData", sep=""))

###########################################################################
