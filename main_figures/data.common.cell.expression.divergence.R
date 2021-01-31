###########################################################################################

options(stringsAsFactors = FALSE)

source("parameters.R") ## pathFinalData are defined based on the user name

path_exp <- paste(pathFinalData, "SupplementaryDataset6/", sep="")
path_evol <- paste(pathFinalData, "SupplementaryDataset7/", sep="")

###########################################################################################

load(paste(pathFigures, "RData/data.gene.annotations.RData", sep=""))
load(paste(pathFigures, "RData/data.enhancer.statistics.RData", sep=""))

ortho <- read.table(paste(path_evol, "human/gene_orthology/human2mouse_orthologue_dNdS.txt", sep="/"), h=T, sep="\t")

exp_human = read.table(paste(path_exp, "human/TPMValues_CommonCellTypes.txt", sep=""), h=T, stringsAsFactors=F, row.names = 1, sep="\t")
exp_mouse = read.table(paste(path_exp, "mouse/TPMValues_CommonCellTypes.txt", sep=""), h=T, stringsAsFactors=F, row.names = 1, sep="\t")

expdiv_cells=read.table(paste(path_exp, "expression_divergence/ExpressionDivergence_CellTypes_MeanTPM.txt",sep=""), h=T, stringsAsFactors=F, sep="\t")
colnames(expdiv_cells) = c("IDHuman", "IDMouse", "adipo", "Bcell", "ESC")
rownames(expdiv_cells) = expdiv_cells$IDHuman
expdiv_cells <- expdiv_cells[complete.cases(expdiv_cells), ] 

cells <- c("Bcell", "ESC", "adipo")
ESC.common = list("human"= c("hESC"), "mouse"=c("ESC", "ESC_18", "ESC_wild"))
adipo.common = list("human"= c("pre_adipo"), "mouse"=c("preadip_D0", "preadip_D2", "preadip_4H"))
Bcell.common = list("human"= c("TB", "NB"), "mouse"=c("preB_aged", "preB_young"))
common.cell <- list("ESC"=ESC.common, "adipo"=adipo.common, "Bcell"=Bcell.common)

enh = "ENCODE"

###########################################################################################
################################### Expression divergence #################################
## select genes: protein-coding genes
annot=gene.annot[["human"]]
coding_genes=annot$GeneID[which(annot$GeneBiotype=="protein_coding")]
genes=intersect(rownames(expdiv_cells), coding_genes)
expdiv_cells = expdiv_cells[genes,]

for (cell in cells){
  ##Expression in mouse
  samples = grep(cell, colnames(exp_human), value=T)
  expdiv_cells[[paste0(cell, "_human_MeanRPKM")]] <- apply(exp_human[expdiv_cells[,"IDHuman"], samples], 1, function(x) mean(x))

  if (cell == "Bcell"){
    cell_name="B"
  } else{
    cell_name = cell # Because Bcell in mouse are intitulate B1, B2...
  }
  
  samples = grep(cell_name, colnames(exp_mouse), value=T)

  expdiv_cells[[paste0(cell, "_mouse_MeanRPKM")]] <- apply(exp_mouse[expdiv_cells[,"IDMouse"], samples], 1, function(x) mean(x))

  global.mean <- (expdiv_cells[[paste0(cell, "_human_MeanRPKM")]] + expdiv_cells[[paste0(cell, "_mouse_MeanRPKM")]])/2
  expdiv_cells[[paste0(cell, "_globalMeanRPKM")]] <- global.mean

  lm1 = lm(1-expdiv_cells[,cell]~log2(expdiv_cells[[paste0(cell, "_globalMeanRPKM")]]+1))

  expdiv_cells[[paste0(cell, "_ExpressionConservation")]] = 1-expdiv_cells[,cell]
  expdiv_cells[[paste0(cell, "_ResidualExpressionConservation")]] = lm1$residuals
}

save(expdiv_cells, file = paste(pathFigures, "/RData/data.common.cells.expdiv.RData", sep=""))

############################################################################################
