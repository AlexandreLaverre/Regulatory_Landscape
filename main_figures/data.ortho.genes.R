###########################################################################

source("parameters.R")

pathEvolution=paste(pathFinalData, "SupplementaryDataset7/", sep="")
pathGenes=paste(pathFinalData, "SupplementaryDataset3/genes/", sep="")

###########################################################################

genes.human=read.table(paste(pathGenes, "human_genes_Ensembl94.txt", sep=""), h=T, stringsAsFactors=F, sep="\t")
genes.mouse=read.table(paste(pathGenes, "mouse_genes_Ensembl94.txt", sep=""), h=T, stringsAsFactors=F, sep="\t")

pc.human=genes.human$GeneID[which(genes.human$GeneBiotype=="protein_coding" & genes.human$Chr%in%c(as.character(1:22), "X", "Y"))]
pc.mouse=genes.mouse$GeneID[which(genes.mouse$GeneBiotype=="protein_coding" & genes.mouse$Chr%in%c(as.character(1:19), "X", "Y"))]

###########################################################################

## only for human and mouse

all.ortho=read.table(paste(pathEvolution, "human/gene_orthology/human2mouse_orthologue_dNdS.txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
colnames(all.ortho)=c("IDHuman", "IDMouse", "HomologyType", "dN", "dS")

all.ortho$dN=as.numeric(all.ortho$dN)
all.ortho$dS=as.numeric(all.ortho$dS)
all.ortho$dNdS=all.ortho$dN/all.ortho$dS

###########################################################################

ortho=all.ortho[which(all.ortho$HomologyType=="ortholog_one2one" & all.ortho$IDHuman%in%pc.human & all.ortho$IDMouse%in%pc.mouse), c("IDHuman", "IDMouse", "dN", "dS", "dNdS")]
colnames(ortho)=c("human", "mouse", "dN", "dS", "dNdS")

###########################################################################

save(ortho, file=paste(pathFigures, "RData/data.ortho.genes.RData", sep=""))

###########################################################################
