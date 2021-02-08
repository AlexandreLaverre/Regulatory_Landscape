library(Hmisc)
options(stringsAsFactors = FALSE)

#setwd("/home/laverre/Data/Regulatory_landscape/scripts/main_figures")
source("parameters.R") ## pathFinalDataFinalDatas are defined based on the user name

path_exp <- paste(pathFinalData, "SupplementaryDataset6/", sep="")
path_evol <- paste(pathFinalData, "SupplementaryDataset7/", sep="")

load(paste(pathFigures, "RData/data.gene.annotations.RData", sep=""))

#######################################################################################################

for (sp in c("human", "mouse")){
  if (sp == "human"){sp_name="Human"}else{sp_name="Mouse"}
  
  expdiv=read.table(paste0(pathFinalData, "SupplementaryDataset6/expression_divergence/ExpressionDivergence_CardosoMoreira2019_SomaticOrgans.txt"), h=T, stringsAsFactors=F, sep="\t")
  rownames(expdiv)=expdiv[[paste0("ID", sp_name)]]
  
  expdiv$classTau=cut2(expdiv[[paste0("Tau", sp_name)]], g=4, include.lowest=T)
  
  #select genes: protein-coding genes
  annot=gene.annot[[sp]]
  coding_genes=annot$GeneID[which(annot$GeneBiotype=="protein_coding")]
  genes=intersect(rownames(expdiv), coding_genes)
  
  expdiv = expdiv[genes,]
  
  samples.human=setdiff(grep("^Human_", colnames(expdiv), value=T), "Human_MeanRPKM")
  samples.mouse=setdiff(grep("^Mouse_", colnames(expdiv), value=T), "Mouse_MeanRPKM")
  shh.human=as.numeric(expdiv["ENSG00000164690", samples.human])
  shh.mouse=as.numeric(expdiv["ENSG00000164690", samples.mouse])
  
  ####################################################################################################
  # Output
  save(expdiv, samples.human, samples.mouse, shh.human, shh.mouse,
       file = paste(pathFigures, "/RData/data.", sp, ".CM2019.SomaticOrgans.expdiv.Rdata", sep=""))
}

