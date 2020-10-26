################################################################################################################################
library(gsubfn)
library(data.table)
library(Hmisc)
options(stringsAsFactors = FALSE)

#setwd("/home/laverre/Data/Regulatory_landscape/scripts/main_figures")
source("parameters.R") ## pathFinalDataFinalDatas are defined based on the user name
path_exp <- paste(pathFinalData, "SupplementaryDataset6/", sep="")
path_evol <- paste(pathFinalData, "SupplementaryDataset7/", sep="")

ref_sp = "human" 
if (ref_sp == "human"){target_sp = "mouse"}else{target_sp = "human"}

cells <- c("Bcell", "ESC", "adipo")
enh = "ENCODE"
dataset.colors=c("firebrick1", "forestgreen", "navy")
names(dataset.colors) = cells

################################################################################################################
############################################## Datas ###########################################################

ortho <- read.table(paste(path_evol, ref_sp, "/gene_orthology/human2mouse_orthologue_dNdS.txt", sep="/"), h=T, sep="\t")
rownames(ortho) <- ortho$GenestableID
ortho$dNdS <- ortho$dN/ortho$dS
ortho <- ortho[which(!is.na(ortho$dNdS) & ortho$dNdS < 50),]

expdiv=read.table(paste(path_exp, "expression_divergence/ExpressionDivergence_CellTypes_MedianTPM.txt",sep=""), h=T, stringsAsFactors=F, sep="\t")
colnames(expdiv) = c("IDHuman", "IDMouse", "adipo", "Bcell", "ESC")
expdiv <- expdiv[complete.cases(expdiv), ] # to remove if we find a solution to NA

exp_human = read.table(paste(path_exp, "human/TPMValues_CommonCellTypes.txt", sep=""), h=T, stringsAsFactors=F, row.names = 1, sep="\t")
exp_mouse = read.table(paste(path_exp, "mouse/TPMValues_CommonCellTypes.txt", sep=""), h=T, stringsAsFactors=F, row.names = 1, sep="\t")

if(ref_sp == "human"){rownames(expdiv)=expdiv$IDHuman}else{rownames(expdiv)=expdiv$IDMouse}

nb_class=20
for (cell in cells){
  samples = grep(cell, colnames(exp_human), value=T)
  expdiv[[paste0("human_", cell, "_Mean")]] <- apply(exp_human[expdiv[,"IDHuman"], samples], 1, function(x) mean(x))
  expdiv[[paste0("human_", cell, "_expClass")]] <- cut2(expdiv[[paste0("human_", cell, "_Mean")]], g=nb_class)
  levels(expdiv[[paste0("human_", cell, "_expClass")]]) <- 1:nb_class
  
  if (cell == "Bcell"){cell_name="B"}else{cell_name = cell} # Because Bcell in mouse are intitulate B1, B2...
  samples = grep(cell_name, colnames(exp_mouse), value=T)
  expdiv[[paste0("mouse_", cell, "_Mean")]] <- apply(exp_mouse[expdiv[,"IDMouse"], samples], 1, function(x) mean(x))
  expdiv[[paste0("mouse_", cell, "_expClass")]] <- cut2(expdiv[[paste0("mouse_", cell, "_Mean")]], g=nb_class)
  levels(expdiv[[paste0("mouse_", cell, "_expClass")]]) <- 1:nb_class
  
  expdiv[[paste0(cell, "_Mean")]] <- (expdiv[[paste0("human_", cell, "_Mean")]] + expdiv[[paste0("mouse_", cell, "_Mean")]])/2
  
  lm1 = lm(expdiv[,cell]~log2(expdiv[[paste0(cell, "_Mean")]]+1))
  expdiv[[paste0(cell, "_ResidualExpressionDivergence")]] = lm1$residuals
  
}

par(mfrow=c(1,3))
for (cell in cells){
  mat <- table(expdiv[[paste0("human_", cell, "_expClass")]], expdiv[[paste0("mouse_", cell, "_expClass")]])
  heatmap.2(mat, trace="none", Rowv = F, Colv=F, scale="none", dendrogram="none", col =colorRampPalette(c("white","navyblue", "black"))(200),
            key.xlab="Gene density", key.title = "",
            main=cell, xlab = "Mouse Expression Quantile", ylab ="Human Expression Quantile")
}

##############################################################################################################################
data_cell <- list()

for (cell in cells){
  regland = read.table(paste(path_evol, ref_sp, "/evolution_summary_by_gene/", enh, "/original_evolution_summary_", cell, "_0.4.txt",sep=""), h=T, stringsAsFactors=F, sep="\t", row.names = 1)
  common=intersect(rownames(expdiv), rownames(regland))
  expdiv_cell=expdiv[common,]
  regland=regland[common,]
  
  # Made decile of nb enhancers
  regland$class_nb_contact <- cut2(regland$nb_total, g=5)
  
  regland$ratio_cons_seq = regland$nb_seq_conserv/regland$nb_total
  regland$ratio_cons_int = ifelse(regland$nb_seq_conserv > 0, regland$nb_contact_conserv/regland$nb_seq_conserv, 0)
  
  regland$ratio_cons_synt = ifelse(regland$nb_seq_conserv > 5 & regland$nb_seq_conserv < 50, regland$nb_synt2M_conserv/regland$nb_seq_conserv, NA)
  regland$ratio_cons_int = ifelse(regland$nb_seq_conserv > 5 & regland$nb_seq_conserv < 50, regland$nb_contact_conserv/regland$nb_seq_conserv, NA)
  
  regland$class_align_score=cut2(regland$med_align_score, g=5, include.lowest=T)
  regland$class_cons_seq=cut(regland$ratio_cons_seq, breaks=c(0, 0.001, 0.25, 0.50, 0.75, 1), include.lowest=T)
  regland$class_cons_int=cut(regland$ratio_cons_int,  breaks=c(0, 0.001, 0.25, 0.50, 0.75, 1), include.lowest=T)
  
  regland$divergence <- expdiv_cell[[paste0(cell)]] # expdiv[[paste0(cell, "_ResidualExpressionDivergence")]]
  regland$human_class <- expdiv_cell[[paste0("human_", cell, "_expClass")]]
  regland$mouse_class <- expdiv_cell[[paste0("mouse_", cell, "_expClass")]]
  
  data_cell[[cell]] <- regland
}

###########################################################################################################################
par(mfrow=c(1,3))
for (cell in cells){
  mat=matrix(rep(NA, nb_class^2), nrow=nb_class)
  mat_verif=matrix(rep(NA, nb_class^2), nrow=nb_class)
  rownames(mat)=1:nb_class
  colnames(mat)=1:nb_class
  rownames(mat_verif)=1:nb_class
  colnames(mat_verif)=1:nb_class
  
  for(i in 1:nb_class){
    class_human=levels(expdiv[[paste0("human_", cell, "_expClass")]])[i]
    for(j in 1:nb_class){
      class_mouse=levels(expdiv[[paste0("mouse_", cell, "_expClass")]])[j]
      
      genes = data_cell[[cell]][which(data_cell[[cell]]$human_class == class_human & data_cell[[cell]]$mouse_class == class_mouse),]
      if (nrow(genes) > 5){
        mat[class_human, class_mouse]=median(genes$ratio_cons_synt)
      }else{mat[class_human, class_mouse]=NA}
      
      mat_verif[class_human, class_mouse]=nrow(genes)
    }
  }
  
  heatmap.2(mat, trace="none", Rowv = F, Colv=F, scale="none", dendrogram="none", 
            col = rev(colorRampPalette(brewer.pal(8, "Spectral"))(25)), na.color="grey95",
            key.xlab="Median Ratio conserved synt", key.title = "",
            main=cell, xlab = "Mouse Expression Quantile", ylab ="Human Expression Quantile")
}

