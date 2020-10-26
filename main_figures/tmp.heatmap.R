################################################################################################################################
library(gsubfn)
library(data.table)
library(Hmisc)
library(gplots)
library(RColorBrewer)
options(stringsAsFactors = FALSE)

#setwd("/home/laverre/Data/Regulatory_landscape/scripts/main_figures")
source("parameters.R") ## pathFinalDataFinalDatas are defined based on the user name
path_exp <- paste(pathFinalData, "SupplementaryDataset6/", sep="")
path_evol <- paste(pathFinalData, "SupplementaryDataset7/", sep="")
load(paste(pathFigures, "RData/data.gene.enhancer.align.RData", sep=""))
load(paste(pathFigures, "RData/data.gene.annotations.RData", sep=""))

ref_sp = "human" 
if (ref_sp == "human"){target_sp = "mouse"}else{target_sp = "human"}

cells <- c("Testis") #c("Bcell", "ESC", "adipo")
enh = "ENCODE"
dataset.colors=c("firebrick1", "forestgreen", "navy")
names(dataset.colors) = cells

nb_class=25

################################################################################################################
############################################## Datas ###########################################################

ortho <- read.table(paste(path_evol, ref_sp, "/gene_orthology/human2mouse_orthologue_dNdS.txt", sep="/"), h=T, sep="\t")
rownames(ortho) <- ortho$GenestableID
ortho$dNdS <- ortho$dN/ortho$dS
ortho <- ortho[which(!is.na(ortho$dNdS) & ortho$dNdS < 50),]

# expdiv=read.table(paste(path_exp, "expression_divergence/ExpressionDivergence_CellTypes_MedianTPM.txt",sep=""), h=T, stringsAsFactors=F, sep="\t")
# colnames(expdiv) = c("IDHuman", "IDMouse", "adipo", "Bcell", "ESC")
# expdiv <- expdiv[complete.cases(expdiv), ] # to remove if we find a solution to NA
# 
# exp_human = read.table(paste(path_exp, "human/TPMValues_CommonCellTypes.txt", sep=""), h=T, stringsAsFactors=F, row.names = 1, sep="\t")
# exp_mouse = read.table(paste(path_exp, "mouse/TPMValues_CommonCellTypes.txt", sep=""), h=T, stringsAsFactors=F, row.names = 1, sep="\t")
# 
# for (cell in cells){
#   samples = grep(cell, colnames(exp_human), value=T)
#   expdiv[[paste0("human_", cell, "_Mean")]] <- apply(exp_human[expdiv[,"IDHuman"], samples], 1, function(x) mean(x))
#   expdiv[[paste0("human_", cell, "_expClass")]] <- cut2(expdiv[[paste0("human_", cell, "_Mean")]], g=nb_class)
#   levels(expdiv[[paste0("human_", cell, "_expClass")]]) <- 1:nb_class
# 
#   if (cell == "Bcell"){cell_name="B"}else{cell_name = cell} # Because Bcell in mouse are intitulate B1, B2...
#   samples = grep(cell_name, colnames(exp_mouse), value=T)
#   expdiv[[paste0("mouse_", cell, "_Mean")]] <- apply(exp_mouse[expdiv[,"IDMouse"], samples], 1, function(x) mean(x))
#   expdiv[[paste0("mouse_", cell, "_expClass")]] <- cut2(expdiv[[paste0("mouse_", cell, "_Mean")]], g=nb_class)
#   levels(expdiv[[paste0("mouse_", cell, "_expClass")]]) <- 1:nb_class
# 
#   expdiv[[paste0(cell, "_Mean")]] <- (expdiv[[paste0("human_", cell, "_Mean")]] + expdiv[[paste0("mouse_", cell, "_Mean")]])/2
# 
#   lm1 = lm(expdiv[,cell]~log2(expdiv[[paste0(cell, "_Mean")]]+1))
#   expdiv[[paste0(cell, "_ResidualExpressionDivergence")]] = lm1$residuals
# }

# Testis_S12 = Testis_P63 for mouse and Testis_youngMidAge for human
expdiv_all=read.table(paste(path_exp, "expression_divergence/ExpressionDivergence_CardosoMoreira2019_OrganStage_MeanTPM.txt",sep=""), h=T, stringsAsFactors=F, sep="\t")
if (ref_sp == "human"){rownames(expdiv_all)=expdiv_all$IDHuman}else{rownames(expdiv_all)=expdiv_all$IDMouse}
expdiv_all <- expdiv_all[complete.cases(expdiv_all), ] # to remove if we find a solution to NA
expdiv <- expdiv_all[,1:2]
expdiv[[cell]] <- expdiv_all$Testis_S12

exp_human = read.table(paste(path_exp, "human/AverageRPKM_CardosoMoreira2019.txt", sep=""), h=T, stringsAsFactors=F, row.names = 1, sep="\t")
exp_mouse = read.table(paste(path_exp, "mouse/AverageRPKM_CardosoMoreira2019.txt", sep=""), h=T, stringsAsFactors=F, row.names = 1, sep="\t")

for (cell in cells){
  expdiv[[paste0("human_", cell, "_Mean")]] <- exp_human[expdiv[,"IDHuman"], paste0(cell, ".youngMidAge")]
  expdiv[[paste0("human_", cell, "_expClass")]] <- cut2(expdiv[[paste0("human_", cell, "_Mean")]], g=nb_class)
  levels(expdiv[[paste0("human_", cell, "_expClass")]]) <- 1:nb_class

  expdiv[[paste0("mouse_", cell, "_Mean")]] <- exp_mouse[expdiv[,"IDMouse"], paste0(cell, ".P63")]
  expdiv[[paste0("mouse_", cell, "_expClass")]] <- cut2(expdiv[[paste0("mouse_", cell, "_Mean")]], g=nb_class)
  levels(expdiv[[paste0("mouse_", cell, "_expClass")]]) <- 1:nb_class

  expdiv[[paste0(cell, "_Mean")]] <- (expdiv[[paste0("human_", cell, "_Mean")]] + expdiv[[paste0("mouse_", cell, "_Mean")]])/2

  lm1 = lm(expdiv[,cell]~log2(expdiv[[paste0(cell, "_Mean")]]+1))
  expdiv[[paste0(cell, "_ResidualExpressionDivergence")]] = lm1$residuals
}


## select genes: protein-coding genes, in the PCHiC data, and in expression stats
annot=gene.annot[[ref_sp]]
coding_genes=annot$GeneID[which(annot$GeneBiotype=="protein_coding")]
coding_genes=intersect(rownames(expdiv), coding_genes)
expdiv = expdiv[coding_genes,]

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
  #regland = read.table(paste(path_evol, ref_sp, "/evolution_summary_by_gene/", enh, "/original_evolution_summary_", cell, "_0.4.txt",sep=""), h=T, stringsAsFactors=F, sep="\t", row.names = 1)
  regland = read.table(paste(path_evol, ref_sp, "/evolution_summary_by_gene/", enh, "/original_evolution_summary_all_0.4.txt",sep=""), h=T, stringsAsFactors=F, sep="\t", row.names = 1)
  
  common=intersect(rownames(expdiv), rownames(regland))
  expdiv_cell=expdiv[common,]
  regland=regland[common,]
  
  # Made decile of nb enhancers
  regland$class_nb_contact <- cut2(regland$nb_total, g=5)
  
  regland$ratio_cons_seq = ifelse(regland$nb_total > 5 & regland$nb_total < 50, regland$nb_seq_conserv/regland$nb_total, NA)
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
# par(mfrow=c(1,3))
# for (cell in cells){
#   mat=matrix(rep(NA, nb_class^2), nrow=nb_class)
#   mat_verif=matrix(rep(NA, nb_class^2), nrow=nb_class)
#   rownames(mat)=1:nb_class
#   colnames(mat)=1:nb_class
#   rownames(mat_verif)=1:nb_class
#   colnames(mat_verif)=1:nb_class
#   enhancers.align = gene.enhancers.align[[ref_sp]][[enh]][[cell]][["real"]]
#   
#   for(i in 1:nb_class){
#     class_human=levels(expdiv[[paste0("human_", cell, "_expClass")]])[i]
#     for(j in 1:nb_class){
#       class_mouse=levels(expdiv[[paste0("mouse_", cell, "_expClass")]])[j]
#       
#       genes = row.names(data_cell[[cell]][which(data_cell[[cell]]$human_class == class_human & data_cell[[cell]]$mouse_class == class_mouse),])
#       if (length(genes) > 5){
#         contacted.enh = unique(enhancers.align[which(enhancers.align$gene %in% genes), c("enhancer", "alignmentscore")])
#         mat[class_human, class_mouse]=mean(contacted.enh$alignmentscore, na.rm=T)
#       }else{mat[class_human, class_mouse]=NA}
#       
#       mat_verif[class_human, class_mouse]=length(genes)
#     }
#   }
#   
#   heatmap.2(mat, trace="none", Rowv = F, Colv=F, scale="none", dendrogram="none", 
#             col = rev(colorRampPalette(brewer.pal(11, "Spectral"))(25)), na.color="grey95",
#             key.xlab="Median Enhancer Alignment Score", key.title = "",
#             main=cell, xlab = "Mouse Expression Quantile", ylab ="Human Expression Quantile")
# }

###########################################################################################################################
###########################################################################################################################
# Other measure
# 
vars = c("divergence", "nb_total", "ratio_cons_seq", "med_align_score", "ratio_cons_int")
names(vars) = c("Median Expression Divergence", "Median Nb Enhancers", "Median Ratio conserved Enh", "Median Alignment Score", "Median Ratio conserved contacts")
for (var in vars){
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
          mat[class_human, class_mouse]=median(genes[[var]], na.rm = T)
        }else{mat[class_human, class_mouse]=NA}
        
        mat_verif[class_human, class_mouse]=nrow(genes)
      }
    }
    
    heatmap.2(mat, trace="none", Rowv = F, Colv=F, scale="none", dendrogram="none",
              col = rev(colorRampPalette(brewer.pal(11, "Spectral"))(25)), na.color="grey95",
              key.xlab=names(which(vars==var)), key.title = "",
              main=cell, xlab = "Mouse Expression Quantile", ylab ="Human Expression Quantile")
  }
}


