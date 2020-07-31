################################################################################################################################
library(gsubfn)
library(data.table)
options(stringsAsFactors = FALSE)

#setwd("/home/laverre/Data/Regulatory_landscape/scripts/main_figures")
source("parameters.R") ## pathFinalDataFinalDatas are defined based on the user name

ref_sp = "human" 
if (ref_sp == "human"){target_sp = "mouse"}else{target_sp = "human"}


cells <- c("Bcell", "ESC", "adipo")
enhancers = c("FANTOM5", "ENCODE")
path_exp <- paste(pathFinalData, "SupplementaryDataset6/", sep="")
path_evol <- paste(pathFinalData, "SupplementaryDataset7/", sep="")

################################################################################################################################
################################################################################################################################
# Datas
ortho <- read.table(paste(pathFinalData, "SupplementaryDataset3/human2mouse_ortholog_one2one_genes_Ensembl94", sep="/"), h=T, sep="\t")
rownames(ortho) <- ortho$GenestableID
ortho$dNdS <- ortho$dN/ortho$dS
ortho <- ortho[which(!is.na(ortho$dNdS) & ortho$dNdS < 50),]

expdiv=read.table(paste(path_exp, "expression_divergence/ExpressionDivergence_CellTypes_MeanTPM.txt",sep=""), h=T, stringsAsFactors=F, sep="\t")
colnames(expdiv) = c("IDHuman", "IDMouse", "adipo", "Bcell", "ESC")
expdiv <- expdiv[complete.cases(expdiv), ] # to remove if we find a solution to NA

exp_human = read.table(paste(path_exp, "human/TPMValues_CommonCellTypes.txt", sep=""), h=T, stringsAsFactors=F, row.names = 1, sep="\t")
exp_mouse = read.table(paste(path_exp, "mouse/TPMValues_CommonCellTypes.txt", sep=""), h=T, stringsAsFactors=F, row.names = 1, sep="\t")

pdf(paste(pathFigures, ref_sp, "_expression_level_vs_expression_divergence.pdf", sep=""), width=8.5, height=9.5)
par(mfrow=c(2,2))
for (cell in cells){
  samples = grep(cell, colnames(exp_human), value=T)
  expdiv[[paste0("human_", cell, "_Mean")]] <- apply(exp_human[expdiv[,"IDHuman"], samples], 1, function(x) mean(x))
  
  if (cell == "Bcell"){cell_name="B"}else{cell_name = cell} # Because Bcell in mouse are intitulate B1, B2...
  samples = grep(cell_name, colnames(exp_mouse), value=T)
  expdiv[[paste0("mouse_", cell, "_Mean")]] <- apply(exp_mouse[expdiv[,"IDMouse"], samples], 1, function(x) mean(x))
  
  expdiv[[paste0(cell, "_Mean")]] <- (expdiv[[paste0("human_", cell, "_Mean")]] + expdiv[[paste0("mouse_", cell, "_Mean")]])/2
  
  ## Correlation between expression divergence and contacts conservation ratio
  expidv_cell <- expdiv #[which(log2(expdiv[[paste0(cell, "_Mean")]]+1)),]
  R=cor(log2(expidv_cell[[paste0(cell, "_Mean")]]+1), expidv_cell[,cell], method="pearson")
  rho=cor(log2(expidv_cell[[paste0(cell, "_Mean")]]+1), expidv_cell[,cell], method="spearman")
  
  smoothScatter(log2(expidv_cell[[paste0(cell, "_Mean")]]+1), expidv_cell[,cell], main=cell, 
                xlab="Expression level log2(mean(tpm))", ylab="Expression Divergence")
  
  mtext(paste("Pearson's R = ", round(R, digits=2), ", Spearman's rho = ",round(rho, digits=2),sep=""), side=3, line=0.5)
  x=log2(expidv_cell[[paste0(cell, "_Mean")]]+1)
  y=expidv_cell[,cell]
  abline(lm(y~x), col="red")
  
  lm1 = lm(expidv_cell[,cell]~log2(expidv_cell[[paste0(cell, "_Mean")]]+1))
  expdiv[[paste0(cell, "_ResidualExpressionDivergence")]] = lm1$residuals
}

dev.off()

if(ref_sp == "human"){rownames(expdiv)=expdiv$IDHuman}else{rownames(expdiv)=expdiv$IDMouse}

################################################################################################################################
################################################################################################################################

load_data <- function(enh, cell){
  regland_all = read.table(paste(path_evol, ref_sp, "/evolution_summary_by_gene/", enh, "/original_evolution_summary_all_sample.txt",sep=""), h=T, stringsAsFactors=F, sep="\t", row.names = 1)
  regland_all = regland_all[which(regland_all$nb_total > 0),]
  
  regland = read.table(paste(path_evol, ref_sp, "/evolution_summary_by_gene/", enh, "/original_evolution_summary_", cell, ".txt",sep=""), h=T, stringsAsFactors=F, sep="\t", row.names = 1)
  regland = regland[which(rownames(regland) %in% rownames(regland_all)),] # Get only gene contacted by at least 1 enh in all samples
  common=intersect(rownames(expdiv), rownames(regland))
  expdiv=expdiv[common,]
  regland=regland[common,]
  
  regland$ratio_cons_seq = regland$nb_seq_conserv/regland$nb_total
  regland$ratio_cons_int = ifelse(regland$nb_seq_conserv > 0, regland$nb_contact_conserv/regland$nb_seq_conserv, 0)
  
  if (enh %in% c("ENCODE", "RoadMap")){breaks=c(1,5,11,25,50, max(regland$nb_total))
  breaks_names = c("1-5","6-10","11-25", "26-50", ">50")}
  
  if (enh == "FANTOM5"){breaks=c(1,2,5,10,15, max(regland$nb_total))
  breaks_names = c("1","2-5", "6-10","11-15", ">15")}
  
  if (enh == "GRO_seq"){breaks=c(1,2,5,10,15, max(regland$nb_total))
  breaks_names = c("1-2","3-5", "6-15","16-25", ">25")}
  
  regland$class_nb_contact=cut(regland$nb_total, breaks=breaks, include.lowest=T)
  regland$class_complexity = ifelse(regland$nb_total < summary(regland$nb_total)[2], "low", "high") # first quantile 
  regland$class_complexity <- factor(regland$class_complexity, levels = c("low", "high"))
  regland$class_cons_seq=cut(regland$ratio_cons_seq, breaks=c(0, 0.001, 0.25, 0.50, 0.75, 1), include.lowest=T)
  regland$class_cons_int=cut(regland$ratio_cons_int,  breaks=c(0, 0.001, 0.25, 0.50, 0.75, 1), include.lowest=T)
  
  return(list(regland, breaks_names))
}


################################################################################################################################
############################################# Gene Level Expression vs Number of enhancers ###################################
pdf(paste(pathFigures, ref_sp, "_nb_enhancer_vs_expression_level.pdf", sep=""), width=8.5, height=9.5)
par(mfrow=c(2,2))
for (cell in cells){
  for (enh in enhancers){
    #pdf(paste(path, "Main_figures/Figure6_human_", enh, "_expression.pdf", sep=""), width=8.5, height=9.5
    datas <- load_data(enh, cell)
    regland <- datas[[1]]
    breaks_names <- datas[[2]]
    expdiv_cell <- expdiv[rownames(regland),]

    boxplot(log2(expdiv_cell[[paste(ref_sp, cell, "Mean", sep="_")]]+1)~regland$class_nb_contact,
            outline=F, notch=T, axes=F, xlab="", ylab="", border="navy", boxwex=0.5, main=paste0(enh, " in ", cell))
    axis(side=2)
    axis(side=1, at=1:length(breaks_names), labels=breaks_names)
    box()
    mtext("Average expression level (log2 TPM)", side=2, line=2.75)
    mtext("Number of contacted enhancers", side=1, line=2.75)

    ## Correlation 
    R=cor(log2(expdiv_cell[[paste(ref_sp, cell, "Mean", sep="_")]]+1), regland$nb_total, method="pearson")
    rho=cor(log2(expdiv_cell[[paste(ref_sp, cell, "Mean", sep="_")]]+1), regland$nb_total, method="spearman")
    
    smoothScatter(log2(expdiv_cell[[paste(ref_sp, cell, "Mean", sep="_")]]+1), regland$nb_total, main=paste(enh, "in", cell),
                  xlab=paste(ref_sp, "gene expression log2(TPM)", sep=" "), ylab="Number of contacted enhancer", cex.lab=1.2)
    
    mtext(paste("Pearson's R = ", round(R, digits=2), ", Spearman's rho = ", round(rho, digits=2),sep=""), side=3, line=0.5, cex=0.8)
    x=log2(expdiv_cell[[paste(ref_sp, cell, "Mean", sep="_")]]+1)
    y=regland$nb_total
    abline(lm(y~x), col="red")
    
  }
}
dev.off()

################################################################################################################################
############################################# Expression divergence vs Number of enhancers  ###################################
pdf(paste(pathFigures, ref_sp, "_nb_enhancer_vs_expression_divergence_boxplot.pdf", sep=""), width=8.5, height=9.5)
par(mfrow=c(2,2))
for (cell in cells){
  for (enh in enhancers){
    datas <- load_data(enh, cell)
    regland <- datas[[1]]
    breaks_names <- datas[[2]]
    expdiv_cell <- expdiv[rownames(regland),]
    
    # Expression Divergence vs Number of enhancers
    boxplot(expdiv_cell[[paste0(cell, "_ResidualExpressionDivergence")]]~regland$class_nb_contact,
            outline=F, notch=T, axes=F, xlab="", ylab="", border="navy", boxwex=0.5, main=paste0(enh, " in ", cell))
    axis(side=2)
    axis(side=1, at=1:length(breaks_names), labels=breaks_names)
    box()
    mtext("Residual expression divergence", side=2, line=2.75)
    mtext("Number of contacted enhancers", side=1, line=2.75)
    
    # Correlation
    R=cor(expdiv_cell[[paste0(cell, "_ResidualExpressionDivergence")]], regland$nb_total, method="pearson")
    rho=cor(expdiv_cell[[paste0(cell, "_ResidualExpressionDivergence")]], regland$nb_total, method="spearman")
    
    smoothScatter(expdiv_cell[[paste0(cell, "_ResidualExpressionDivergence")]], regland$nb_total, main=paste(enh, "in", cell),
                  xlab=paste(ref_sp, "gene expression log2(TPM)", sep=" "), ylab="Number of contacted enhancer", cex.lab=1.2)
    
    mtext(paste("Pearson's R = ", round(R, digits=2), ", Spearman's rho = ", round(rho, digits=2),sep=""), side=3, line=0.5, cex=0.8)
    x=expdiv_cell[[paste0(cell, "_ResidualExpressionDivergence")]]
    y=regland$nb_total
    abline(lm(y~x), col="red")
  }
}
dev.off()

################################################################################################################################
############################################# Expression divergence vs Number of enhancers  ###################################
for (enh in enhancers){
  pdf(paste(pathFigures, ref_sp, "_", enh, "_regulatory_divergence_vs_expression_divergence.pdf", sep=""), width=8.5, height=9.5)
  par(mfrow=c(2,2))
  for (cell in cells){
    regland <- load_data(enh, cell)[[1]]
    
    # Genes with at least some enhancers
    if (enh == "FANTOM5"){min_enh = 2}else{min_enh = 5}
    regland = regland[which(regland$nb_total > min_enh),]

    expdiv_cell <- expdiv
    common=intersect(rownames(expdiv_cell), rownames(regland))
    expdiv_cell=expdiv_cell[common,]
    regland=regland[common,]
    
    # Expression Divergence vs Number of conserved enhancers
    breaks = c(">0.1", "0.1-25", "26-50", "51-75", ">75")

    boxplot(expdiv_cell[[paste0(cell)]]~regland$class_cons_seq, main=paste0(enh, " in ", cell),
            xaxt = "n", yaxt='n',xlab="", ylab="", outline=F, notch=T, boxwex=0.6, outcex=0.2)
    
    axis(side=2)
    axis(side=1, at=1:length(breaks), labels=breaks)
    mtext("Expression divergence", side=2, line=2.75)
    mtext("Enhancers conserved in sequence (%)", side=1, line=2.75)
    
    # Expression Divergence vs Number of conserved contacted enhancers
    boxplot(expdiv_cell[[paste0(cell)]]~regland$class_cons_int, 
            xaxt = "n", yaxt='n', xlab="", ylab="", outline=F, notch=T, boxwex=0.6, outcex=0.2)
    
    axis(side=2)
    axis(side=1, at=1:length(breaks_names), labels=breaks)
    mtext("Residual expression divergence", side=2, line=2.75)
    mtext("Enhancers conserved in contact (%)", side=1, line=2.75)

    ## Correlation between expression divergence and sequence conservation ratio
    R=cor(expdiv_cell[[paste0(cell)]], regland$ratio_cons_seq, method="pearson")
    rho=cor(expdiv_cell[[paste0(cell)]], regland$ratio_cons_seq, method="spearman")
    
    smoothScatter(expdiv_cell[[paste0(cell)]], regland$ratio_cons_seq, main="", 
                  xlab="ExpressionDivergence", ylab="Enhancers conserved in sequence (%)")
    
    mtext(paste("Pearson's R = ", round(R, digits=2), ", Spearman's rho = ",round(rho, digits=2),sep=""), side=3, line=0.5)
    x=expdiv_cell[[paste0(cell)]]
    y=regland$ratio_cons_seq
    abline(lm(y~x), col="red")
    
    ## Correlation between expression divergence and contacts conservation ratio
    R=cor(expdiv_cell[[paste0(cell)]], regland$ratio_cons_int, method="pearson")
    rho=cor(expdiv_cell[[paste0(cell)]], regland$ratio_cons_int, method="spearman")
    
    smoothScatter(expdiv_cell[[paste0(cell)]], regland$ratio_cons_int, main="", 
                  xlab="ExpressionDivergence", ylab="Enhancers conserved in contact (%)")
    
    mtext(paste("Pearson's R = ", round(R, digits=2), ", Spearman's rho = ",round(rho, digits=2),sep=""), side=3, line=0.5)
    x=expdiv_cell[[paste0(cell)]]
    y=regland$ratio_cons_int
    abline(lm(y~x), col="red")
  }
  dev.off()
  }

################################################################################################################################
############################################# Distribution of top expressed Gene divergence  ###################################
gene <- list()
enh_evol <- list()

for (enh in enhancers){
  enhancers_alignment = read.table(paste(path_evol, ref_sp, "sequence_conservation/enhancers", enh, "Alignments_stats_all_species_nonexonic_ungapped.txt", sep="/"), h=T)
  enhancers_contact = read.table(paste(pathFinalData, "SupplementaryDataset4", ref_sp, enh,  "statistics_contacted_enhancers_original.txt", sep="/"), h=T, stringsAsFactors = F)
  enhancers_contact$enh <-  do.call(paste,c(enhancers_contact[c("chr","start","end")],sep=":"))
  
  for (cell in cells){
    if (cell == "Bcell"){cell_name="Bcell"}
    if (cell == "adipo"){cell_name="pre_adipo"}
    if (cell == "ESC"){cell_name="hESC"}
    
    # Enhancers contacted by gene in cell types
    enhancers_contact_in_cell <- enhancers_contact[which(enhancers_contact[[paste0(cell_name)]] > 0),]
    enh_evol[[cell]] <- enhancers_alignment[which(enhancers_alignment$enh %in% enhancers_contact_in_cell$enh),]$mouse
    
    # More expressed genes 
    thirdquantile <- summary(log2(expdiv[[paste(ref_sp, cell, "Mean", sep="_")]]))[5]
    expdiv_top <- expdiv[which(log2(expdiv[[paste(ref_sp, cell, "Mean", sep="_")]]) > thirdquantile ),]
    gene[[cell]] <- ortho[which(ortho$GenestableID %in% expdiv_top$IDHuman),]$dNdS
    
  }
  
}

pdf(paste(pathFigures, ref_sp, "_genes_vs_contacted_enhancers_in_cell_types.pdf", sep=""), width=12, height=6)
par(mfrow=c(1,3))
boxplot(c(gene[1], gene[2], gene[3]), names=c("Bcell", "ESC", "adipo"), main="More expressed genes (3rd Qu.)",
        ylab="Human to Mouse Alignment Score", notch=T, outline=F)

boxplot(c(enh_evol[1], enh_evol[2], enh_evol[3]), names=c("Bcell", "ESC", "adipo"), main="FANTOM5 contacted",
        ylab="Human to Mouse Alignment Score", notch=T)

boxplot(c(enh_evol[1], enh_evol[2], enh_evol[3]), names=c("Bcell", "ESC", "adipo"), main="ENCODE contacted",
        ylab="Human to Mouse Alignment Score", notch=T)

dev.off()
