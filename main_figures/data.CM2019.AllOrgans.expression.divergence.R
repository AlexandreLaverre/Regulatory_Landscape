library(Hmisc)
options(stringsAsFactors = FALSE)

source("parameters.R") ## pathFinalDataFinalDatas are defined based on the user name

path_exp <- paste(pathFinalData, "SupplementaryDataset6/", sep="")
path_evol <- paste(pathFinalData, "SupplementaryDataset7/", sep="")

load(paste(pathFigures, "RData/data.gene.annotations.RData", sep=""))

#######################################################################################################

for (sp in c("human", "mouse")){
  if (sp == "human"){
    sp_name="Human"
  } else{
    sp_name="Mouse"
  }
  
  expdiv=read.table(paste0(pathFinalData, "SupplementaryDataset6/expression_divergence/ExpressionDivergence_CardosoMoreira2019_AllOrgans.txt"), h=T, stringsAsFactors=F, sep="\t")
  rownames(expdiv)=expdiv[[paste0("ID", sp_name)]]
  
  ##select genes: protein-coding genes
  annot=gene.annot[[sp]]
  coding_genes=annot$GeneID[which(annot$GeneBiotype=="protein_coding")]
  genes=intersect(rownames(expdiv), coding_genes)
  
  expdiv = expdiv[genes,]
  
  samples.human=setdiff(grep("^Human_", colnames(expdiv), value=T), "Human_MeanRPKM")
  samples.mouse=setdiff(grep("^Mouse_", colnames(expdiv), value=T), "Mouse_MeanRPKM")
 
  ####################################################################################################
  # Output
  save(expdiv, samples.human, samples.mouse, file = paste(pathFigures, "/RData/data.", sp, ".CM2019.AllOrgans.expdiv.RData", sep=""))
}

####################################################################################################
