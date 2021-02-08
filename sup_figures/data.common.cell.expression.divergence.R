###########################################################################################

options(stringsAsFactors = FALSE)

source("../main_figures/parameters.R") ## pathFinalData are defined based on the user name

path_exp <- paste(pathFinalData, "SupplementaryDataset6/", sep="")
path_evol <- paste(pathFinalData, "SupplementaryDataset7/", sep="")

###########################################################################################

load(paste(pathFigures, "RData/data.enhancer.statistics.RData", sep=""))
load(paste(pathFigures, "RData/data.ortho.genes.RData", sep=""))
load(paste(pathFigures, "RData/data.sample.info.RData", sep=""))

###########################################################################################

exp_human = read.table(paste(path_exp, "human/TPMValues_CommonCellTypes.txt", sep=""), h=T, stringsAsFactors=F, row.names = 1, sep="\t")
exp_mouse = read.table(paste(path_exp, "mouse/TPMValues_CommonCellTypes.txt", sep=""), h=T, stringsAsFactors=F, row.names = 1, sep="\t")

expdiv_cells=read.table(paste(path_exp, "expression_divergence/ExpressionDivergence_CellTypes_MeanTPM.txt",sep=""), h=T, stringsAsFactors=F, sep="\t")
colnames(expdiv_cells) = c("IDHuman", "IDMouse", "adipo", "Bcell", "ESC")
rownames(expdiv_cells) = expdiv_cells$IDHuman

## previously filtered ortho genes

expdiv_cells=expdiv_cells[which(expdiv_cells$IDHuman%in%ortho[,"human"] & expdiv_cells$IDMouse%in%ortho[,"mouse"]),]

###########################################################################################

## corrected expression divergence

for (cell in c("Bcell", "adipo", "ESC")){
  ## expression in human
  samples.human = grep(cell, colnames(exp_human), value=T)
  expdiv_cells[[paste0(cell, "_human_MeanRPKM")]] <- apply(exp_human[expdiv_cells[,"IDHuman"], samples.human], 1, mean, na.rm=T)

  ## expression in mouse
  if (cell == "Bcell"){
    cell_name="B"
  } else{
    cell_name = cell # Because Bcell in mouse are named B1, B2...
  }
  
  samples.mouse = grep(cell_name, colnames(exp_mouse), value=T)
  expdiv_cells[[paste0(cell, "_mouse_MeanRPKM")]] <- apply(exp_mouse[expdiv_cells[,"IDMouse"], samples.mouse], 1, mean, na.rm=T)

  ## average expression level

  global.mean <- (expdiv_cells[, paste0(cell, "_human_MeanRPKM")] + expdiv_cells[, paste0(cell, "_mouse_MeanRPKM")])/2
  expdiv_cells[, paste0(cell, "_globalMeanRPKM")] <- global.mean

  ## expression conservation 
  expdiv_cells[, paste0(cell, "_ExpressionConservation")] = 1-expdiv_cells[,cell]

  ## correct for mean expression levels - only non-NA values
  
  this.expdiv=expdiv_cells[which(expdiv_cells[,paste0(cell, "_globalMeanRPKM")]>0),]

  lm1 = lm(1-this.expdiv[,cell]~log2(this.expdiv[,paste0(cell, "_globalMeanRPKM")]+1))
  corrected=lm1$residuals
  names(corrected)=rownames(this.expdiv)

  expdiv_cells[, paste0(cell, "_ResidualExpressionConservation")] = rep(NA, dim(expdiv_cells)[1])
  expdiv_cells[names(corrected), paste0(cell, "_ResidualExpressionConservation")]=corrected
  
}

############################################################################################

save(expdiv_cells, file = paste(pathFigures, "/RData/data.common.cells.expdiv.RData", sep=""))

############################################################################################
