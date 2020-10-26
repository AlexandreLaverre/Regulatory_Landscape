options(stringsAsFactors = FALSE)

#setwd("/home/laverre/Data/Regulatory_landscape/scripts/main_figures")
source("parameters.R") ## pathFinalDataFinalDatas are defined based on the user name

path_exp <- paste(pathFinalData, "SupplementaryDataset6/", sep="")
path_evol <- paste(pathFinalData, "SupplementaryDataset7/", sep="")

load(paste(pathFigures, "RData/data.gene.annotations.RData", sep=""))

sp="human"
if (sp == "human"){sp_name="Human"}else{sp_name="Mouse"}

#######################################################################################################
################################### Expression divergence #############################################
expdiv_all=read.table(paste0(pathFinalData, "SupplementaryDataset6/expression_divergence/v1/ExpressionDivergence_CardosoMoreira2019_correlations.txt"), h=T, stringsAsFactors=F, sep="\t")
rownames(expdiv_all)=expdiv_all[[paste0("ID", sp_name)]]

expdiv_all$classTau=cut2(expdiv_all[[paste0("Tau", sp_name)]], g=4, include.lowest=T)
expdiv_all$EuclidianSimilarity = 1-expdiv_all$ResidualExpressionDivergence

#select genes: protein-coding genes
annot=gene.annot[[sp]]
coding_genes=annot$GeneID[which(annot$GeneBiotype=="protein_coding")]
genes=intersect(rownames(expdiv_all), coding_genes)

for (enh in enhancer.datasets[[sp]]){
  lm1 = lm(expdiv_all$CorrelationSpearman~expdiv_all$TauHuman)
  expdiv_all$SpearmanResidual = lm1$residuals
}

expdiv_all = expdiv_all[genes,]

#################################################################################################################################
############################################# Regulatory Divergence #############################################################
regland = list()
expdiv = list()
for (enh in enhancer.datasets[[sp]]){
  file = paste(pathFinalData, "SupplementaryDataset7/", sp, "/evolution_summary_by_gene/", enh, "/original_evolution_summary_all_0.4.txt", sep="")
  regland_enh = read.table(file, h=T, stringsAsFactors=F, sep="\t", row.names = 1)
  
  if (enh == "FANTOM5"){nb_min=2}else{nb_min=5}
  
  regland_enh <- regland_enh[which(regland_enh$nb_total >= nb_min & regland_enh$nb_total <= 100),] 
  common=intersect(rownames(expdiv_all), rownames(regland_enh))
  regland_enh = regland_enh[common,]
  expdiv_enh = expdiv_all[common,]
  
  regland_enh$ratio_cons_seq = regland_enh$nb_seq_conserv/regland_enh$nb_total
  regland_enh$ratio_cons_synt = ifelse(regland_enh$nb_seq_conserv > nb_min & regland_enh$nb_seq_conserv < 100, regland_enh$nb_synt2M_conserv/regland_enh$nb_seq_conserv, NA)
  regland_enh$ratio_cons_int = ifelse(regland_enh$nb_seq_conserv > nb_min & regland_enh$nb_seq_conserv < 100, regland_enh$nb_contact_conserv/regland_enh$nb_seq_conserv, NA)
  
  regland_enh$class_nb_contact=cut2(regland_enh$nb_total, g=5, include.lowest=T)
  regland_enh$class_cons_seq=cut(regland_enh$ratio_cons_seq, breaks=c(0, 0.10, 0.25, 0.5, 0.75, 1), include.lowest=T)
  regland_enh$class_cons_synt=cut(regland_enh$ratio_cons_synt, breaks=c(0, 0.75, 0.99, 1), include.lowest=T)
  
  regland_enh$class_cons_int=cut(regland_enh$ratio_cons_int,  breaks=c(0, 0.01, 0.25, 0.5, 0.75, 1), include.lowest=T)
  regland_enh$class_align_score=cut2(regland_enh$med_align_score, g=5, include.lowest=T)
  
  regland[[enh]] = regland_enh
  expdiv[[enh]] = expdiv_enh
}

samples.human=setdiff(grep("^Human_", colnames(expdiv_all), value=T), "Human_MeanRPKM")
samples.mouse=setdiff(grep("^Mouse_", colnames(expdiv_all), value=T), "Mouse_MeanRPKM")
shh.human=as.numeric(expdiv_all["ENSG00000164690", samples.human])
shh.mouse=as.numeric(expdiv_all["ENSG00000164690", samples.mouse])

#########################################################################################################################
# Output

save(regland, expdiv_all, expdiv, samples.human, samples.mouse, shh.human, shh.mouse,
     file = paste(pathFigures, "/RData/Fig6_", sp, "_all_samples_CardosoMoreira.Rdata", sep=""))
