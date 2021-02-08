######################################################################################################################
library(Hmisc)
setwd("/home/laverre/Data/Regulatory_landscape/scripts/main_figures")

objects=ls()

if(!"pathScripts"%in%objects){
  load=T
  prepare=T
  
  source("parameters.R")
}

##############################################################################
regland_break <- list()
expdiv_reg <- list()
for (sp in c("human", "mouse")){
  
  load(paste(pathFigures, "RData/data.gene.annotations.RData", sep=""))
  load(paste(pathFigures, "RData/data.", sp, ".gene.regland.conservation.RData", sep=""))
  
  if (sp == "human"){sp_name="Human"; sp_target="Mouse"}else{sp_name="Mouse"; sp_target="Human"}
  expdiv=read.table(paste0(pathFinalData, "SupplementaryDataset6/expression_divergence/ExpressionDivergence_CardosoMoreira2019_SomaticOrgans.txt"), h=T, stringsAsFactors=F, sep="\t")
  rownames(expdiv)=expdiv[[paste0("ID", sp_name)]]
  
  #select genes: protein-coding genes
  annot=gene.annot[[sp]]
  coding_genes=annot$GeneID[which(annot$GeneBiotype=="protein_coding")]
  genes=intersect(rownames(expdiv), coding_genes)
  expdiv = expdiv[genes,]
  
  ################################################################################################################################
  enh="ENCODE"
  regland = genes.conservation[[enh]][["obs"]][["all"]]
  regland <- regland[which(regland$seq_conserv >= 5 & regland$seq_conserv <= 100),]
  
  genes=intersect(rownames(expdiv), rownames(regland))
  expdiv_reg[[sp]] = expdiv[genes,]
  
  regland$rearranged_ratio = (regland$seq_conserv-regland$synt_conserv)/regland$seq_conserv
  
  regland_break[[sp]] <- regland[which(regland$rearranged >= 0.75), c("seq_conserv", "synt_conserv", "rearranged_ratio")]
  genes = intersect(rownames(regland_break[[sp]]), rownames(expdiv))
  regland_break[[sp]] = regland_break[[sp]][genes,]
  
  expdiv_break = expdiv[genes,]
  regland_break[[sp]]$ortho_genes <- expdiv_break[[paste("ID", sp_target, sep="")]]
  
  write.table(regland_break[[sp]], file=paste(pathFigures, "/test.", sp, ".genes.rearranged.txt",sep=""), row.names=T, col.names=F, quote=F, sep="\t")
  
  genes = intersect(rownames(regland), rownames(expdiv))
  regland_all = regland[genes,]
  
  write.table(regland_all, file=paste(pathFigures, "/test.", sp, ".genes.all.txt",sep=""), row.names=T, col.names=F, quote=F, sep="\t")
  
}
  

common_genes = intersect(rownames(regland_break[["human"]]), rownames(expdiv_reg[["human"]]))
actual_exp = expdiv_reg[["human"]]
actual_exp$RPKM <- log2(actual_exp$Human_MeanRPKM+1)

par(mfrow=c(2,2))
boxplot(actual_exp[common_genes, "RPKM"], actual_exp[which(rownames(actual_exp) %nin% common_genes), "RPKM"],
        notch=T, outline=F, names=c("Highly rearranged", "Other"), ylab="RPKM", cex.lab=1.2)

boxplot(actual_exp[common_genes, "TauHuman"], actual_exp[which(rownames(actual_exp) %nin% common_genes), "TauHuman"],
        notch=T, names=c("Highly rearranged", "Other"), ylab="Specificity", cex.lab=1.2)

boxplot(actual_exp[common_genes, "CorrelationSpearman"], actual_exp[which(rownames(actual_exp) %nin% common_genes), "CorrelationSpearman"],
        notch=T, outline=F, names=c("Highly rearranged", "Other"), ylab="CorrelationSpearman", cex.lab=1.2)




