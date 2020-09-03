################################################################################################################################################
library(data.table)
options(stringsAsFactors = FALSE)

#setwd("/home/laverre/Data/Regulatory_landscape/scripts/main_figures")
source("parameters.R") ## pathFinalDataFinalDatas are defined based on the user name
pathFigures <- "/home/laverre/Documents/Test_Figure7/"

ref_sp = "human" 
if (ref_sp == "human"){target_sp = "mouse"}else{target_sp = "human"}

# Choose genes within : all ; dvpt ; other
selected_genes = "dvpt"

cells <- c("Bcell", "ESC", "adipo")
enhancers = c("FANTOM5", "ENCODE")

path_exp <- paste(pathFinalData, "SupplementaryDataset6/", sep="")
path_evol <- paste(pathFinalData, "SupplementaryDataset7/", sep="")

########################################################  Datas ########################################################  

exp_ref = read.table(paste(path_exp, ref_sp, "/TPMValues_CommonCellTypes.txt", sep=""), h=T, stringsAsFactors=F, row.names = 1, sep="\t")
exp_target = read.table(paste(path_exp, target_sp, "/TPMValues_CommonCellTypes.txt", sep=""), h=T, stringsAsFactors=F, row.names = 1, sep="\t")
ortho = read.table(paste(pathFinalData, "SupplementaryDataset3/human2mouse_ortholog_one2one_genes_Ensembl94.txt", sep=""), h=T, stringsAsFactors=F, sep="\t")
rownames(ortho) <- ortho[[paste0("ID.", ref_sp)]]

ortho = ortho[which(ortho[[paste0("ID.", ref_sp)]] %in% rownames(exp_ref) & ortho[[paste0("ID.", target_sp)]] %in% rownames(exp_target)),]

############################ A - Correlation between gene expression and regulatory landscape complexity ############################ 
pdf(paste(pathFigures, ref_sp, "_correlation_gene_expression_complexity_in_same_cell.pdf", sep=""))
par(mfrow=c(3,2))

for (cell in cells){
  # Calculate mean TPM for each cell type
  if (ref_sp == "mouse" & cell == "Bcell"){cell_name="B"}else{cell_name = cell} # Because Bcell in mouse are intutilate B1, B2...
  samples = grep(cell_name, colnames(exp_ref), value=T)
  ortho[[paste0(ref_sp, "_", cell)]] <- apply(exp_ref[ortho[[paste0("ID.", ref_sp)]], samples], 1, function(x) mean(x))
  
  for (enh in enhancers){
    regland_sp = read.table(paste(path_evol, ref_sp, "/evolution_summary_by_gene/", enh, "/original_evolution_summary_", cell, ".txt",sep=""), h=T, stringsAsFactors=F, sep="\t", row.names = 1)
    
    ortho_regland = ortho[which(ortho[[paste0("ID.", ref_sp)]] %in% rownames(regland_sp)),]
    ortho_regland$sp_complex <- regland_sp[ortho_regland[[paste0("ID.", ref_sp)]], "nb_total"]
    
    
    ## Correlation 
    R=cor(log2(ortho_regland[[paste0(ref_sp, "_", cell)]]+1), ortho_regland$sp_complex, method="pearson")
    rho=cor(log2(ortho_regland[[paste0(ref_sp, "_", cell)]]+1), ortho_regland$sp_complex, method="spearman")
    
    smoothScatter(log2(ortho_regland[[paste0(ref_sp, "_", cell)]]+1), ortho_regland$sp_complex, main=paste(enh, "in", cell),
                  xlab=paste(ref_sp, "gene expression log2(TPM)", sep=" "), ylab="Number of contacted enhancer", cex.lab=1.2)
    
    mtext(paste("Pearson's R = ", round(R, digits=2), ", Spearman's rho = ", round(rho, digits=2),sep=""), side=3, line=0.5, cex=0.8)
    x=log2(ortho_regland[[paste0(ref_sp, "_", cell)]]+1)
    y=ortho_regland$sp_complex
    abline(lm(y~x), col="red")
  }
}
dev.off()


############################ Correlation between human and mouse gene expression in same cell ############################
pdf(paste(pathFigures, "Human_vs_Mouse_gene_expression_correlation_in_same_cell.pdf", sep=""))
par(mfrow=c(2,2))
for (cell in cells){
  # Calculate mean TPM for each cell type
  if (ref_sp == "mouse" & cell == "Bcell"){cell_name="B"}else{cell_name = cell} # Because Bcell in mouse are intutilate B1, B2...
  samples = grep(cell_name, colnames(exp_ref), value=T)
  ortho[[paste0(ref_sp, "_", cell)]] <- apply(exp_ref[ortho[[paste0("ID.", ref_sp)]], samples], 1, function(x) mean(x))
  
  if (target_sp == "mouse" & cell == "Bcell"){cell_name="B"}else{cell_name = cell} # Because Bcell in mouse are intutilate B1, B2...
  samples = grep(cell_name, colnames(exp_target), value=T)
  ortho[[paste0(target_sp, "_", cell)]] <- apply(exp_target[ortho[[paste0("ID.", target_sp)]],samples], 1, function(x) mean(x))
  
  ## Correlation 
  R=cor(log2(ortho[[paste0(ref_sp, "_", cell)]]+1), log2(ortho[[paste0(target_sp, "_", cell)]]+1), method="pearson")
  rho=cor(log2(ortho[[paste0(ref_sp, "_", cell)]]+1), log2(ortho[[paste0(target_sp, "_", cell)]]+1), method="spearman")
  
  smoothScatter(log2(ortho[[paste0(ref_sp, "_", cell)]]+1), log2(ortho[[paste0(target_sp, "_", cell)]]+1), main=paste(cell),
                xlab=paste(ref_sp, "gene expression log2(TPM)", sep=" "), ylab=paste(target_sp, "gene expression log2(TPM)", sep=" "))
  
  mtext(paste("Pearson's R = ", round(R, digits=2), ", Spearman's rho = ", round(rho, digits=2),sep=""), side=3, line=0.5, cex=0.8)
  x=log2(ortho[[paste0(ref_sp, "_", cell)]]+1)
  y=log2(ortho[[paste0(target_sp, "_", cell)]]+1)
  abline(lm(y~x), col="red")
  
}
dev.off()

############################ Correlation between regulatory landscape complexity ############################ 
pdf(paste(pathFigures, "Regulatory_landscape_complexity_correlation_in_same_cell.pdf", sep=""))
par(mfrow=c(3,2))
for (cell in cells){
  for (enh in enhancers){
    regland_ref_sp = read.table(paste(path_evol, ref_sp, "/evolution_summary_by_gene/", enh, "/original_evolution_summary_", cell, ".txt",sep=""), h=T, stringsAsFactors=F, sep="\t", row.names = 1)
    regland_target_sp = read.table(paste(path_evol, target_sp, "/evolution_summary_by_gene/", enh, "/original_evolution_summary_", cell, ".txt",sep=""), h=T, stringsAsFactors=F, sep="\t", row.names = 1)
    
    # Calcul zscore for ortholog genes
    regland_ref_sp$zscore_complex <- (regland_ref_sp$nb_total-mean(regland_ref_sp$nb_total))/sd(regland_ref_sp$nb_total)
    regland_target_sp$zscore_complex <- (regland_target_sp$nb_total-mean(regland_target_sp$nb_total))/sd(regland_target_sp$nb_total)
    
    ortho_regland = ortho[which(ortho[[paste0("ID.", ref_sp)]] %in% rownames(regland_ref_sp) & ortho[[paste0("ID.", target_sp)]] %in% rownames(regland_target_sp)),]
    ortho_regland$ref_sp_complex <- regland_ref_sp[ortho_regland[[paste0("ID.", ref_sp)]], "zscore_complex"]
    ortho_regland$target_sp_complex <- regland_target_sp[ortho_regland[[paste0("ID.", target_sp)]], "zscore_complex"]
    
    
    ## Correlation 
    R=cor(ortho_regland$ref_sp_complex, ortho_regland$target_sp_complex, method="pearson")
    rho=cor(ortho_regland$ref_sp_complex, ortho_regland$target_sp_complex, method="spearman")
    
    smoothScatter(ortho_regland$ref_sp_complex, ortho_regland$target_sp_complex, main=paste(enh, "in", cell),
                  xlab=paste(ref_sp, "complexity (zscore)", sep=" "), ylab=paste(target_sp, "complexity (zscore)", sep=" "), cex.lab=1.2)
    
    mtext(paste("Pearson's R = ", round(R, digits=2), ", Spearman's rho = ", round(rho, digits=2),sep=""), side=3, line=0.5, cex=0.8)
    x=ortho_regland$ref_sp_complex
    y=ortho_regland$target_sp_complex
    abline(lm(y~x), col="red")
  }
}
dev.off()


regland_sp <- regland_sp[which(regland_sp$nb_total > 10),]
regland_sp$ratio_cons_seq = regland_sp$nb_seq_conserv/regland_sp$nb_total

regland_sp <- regland_sp[which(regland_sp$nb_seq_conserv > 10),]
regland_sp$ratio_cons_int = regland_sp$nb_contact_conserv/regland_sp$nb_seq_conserv


R=cor(regland_sp$nb_seq_conserv, regland_sp$ratio_cons_int, method="pearson")
rho=cor(regland_sp$nb_seq_conserv, regland_sp$ratio_cons_int, method="spearman")

smoothScatter(regland_sp$nb_seq_conserv, regland_sp$ratio_cons_int, main=paste(enh, "in", cell),
              xlab="Nb seq", ylab="Ratio conserv", cex.lab=1.2)

mtext(paste("Pearson's R = ", round(R, digits=2), ", Spearman's rho = ", round(rho, digits=2),sep=""), side=3, line=0.5, cex=0.8)
x=regland_sp$nb_seq_conserv
y=regland_sp$ratio_cons_int
abline(lm(y~x), col="red")














############################ Correlation between gene expression and ratio of conserved enhancers ############################ 
pdf(paste(pathFigures, ref_sp, "_correlation_gene_expression_conserved_enhancers_in_same_cell.pdf", sep=""))
par(mfrow=c(3,2))
for (cell in cells){
  # Calculate mean TPM for each cell type
  if (ref_sp == "mouse" & cell == "Bcell"){cell_name="B"}else{cell_name = cell} # Because Bcell in mouse are intitulate B1, B2...
  samples = grep(cell_name, colnames(exp_ref), value=T)
  ortho[[paste0(ref_sp, "_", cell)]] <- apply(exp_ref[ortho[[paste0("ID.", ref_sp)]],samples], 1, function(x) mean(x))
  
  for (enh in enhancers){
    regland_sp = read.table(paste(path_evol, ref_sp, "/evolution_summary_by_gene/", enh, "/original_evolution_summary_", cell, ".txt",sep=""), h=T, stringsAsFactors=F, sep="\t", row.names = 1)
    regland_sp = regland_sp[which(regland_sp$nb_total > 0),]
    regland_sp$ratio_cons_seq = regland_sp$nb_seq_conserv*100/regland_sp$nb_total
    
    ortho_regland = ortho[which(ortho[[paste0("ID.", ref_sp)]] %in% rownames(regland_sp)),]
    ortho_regland$sp_seq_conserv <- regland_sp[ortho_regland[[paste0("ID.", ref_sp)]], "ratio_cons_seq"]
    
    ## Correlation 
    R=cor(log2(ortho_regland[[paste0(ref_sp, "_", cell)]]+1), ortho_regland$sp_seq_conserv, method="pearson")
    rho=cor(log2(ortho_regland[[paste0(ref_sp, "_", cell)]]+1), ortho_regland$sp_seq_conserv, method="spearman")
    
    smoothScatter(log2(ortho_regland[[paste0(ref_sp, "_", cell)]]+1), ortho_regland$sp_seq_conserv, main=paste(enh, "in", cell),
                  xlab=paste(ref_sp, "gene expression log2(TPM)", sep=" "), ylab="Enhancer conserved in sequence (%)", cex.lab=1.2)
    
    mtext(paste("Pearson's R = ", round(R, digits=2), ", Spearman's rho = ", round(rho, digits=2),sep=""), side=3, line=0.5, cex=0.8)
    x=log2(ortho_regland[[paste0(ref_sp, "_", cell)]]+1)
    y=ortho_regland$sp_seq_conserv
    abline(lm(y~x), col="red")
  }
}
dev.off()

############################ Gene expression vs enhancers conserved in contacts ############################ 
pdf(paste(pathFigures, ref_sp, "_expression_level_vs_enhancers_conserved_in_contacts.pdf", sep=""))
par(mfrow=c(3,2))

for (cell in cells){
  # Calculate mean TPM for each cell type
  if (ref_sp == "mouse" & cell == "Bcell"){cell_name="B"}else{cell_name = cell} # Because Bcell in mouse are intutilate B1, B2...
  samples = grep(cell_name, colnames(exp_ref), value=T)
  ortho[[paste0(ref_sp, "_", cell)]] <- apply(exp_ref[ortho[[paste0("ID.", ref_sp)]],samples], 1, function(x) mean(x))
  
  for (enh in enhancers){
    regland_sp = read.table(paste(path_evol, ref_sp, "/evolution_summary_by_gene/", enh, "/original_evolution_summary_", cell, ".txt",sep=""), h=T, stringsAsFactors=F, sep="\t", row.names = 1)
    regland_sp = regland_sp[which(regland_sp$nb_total > 0),]
    regland_sp$ratio_cons_int =  ifelse(regland_sp$nb_seq_conserv > 0, regland_sp$nb_contact_conserv*100/regland_sp$nb_seq_conserv, 0)
    
    ortho_regland = ortho[which(ortho[[paste0("ID.", ref_sp)]] %in% rownames(regland_sp)),]
    ortho_regland$sp_contact_conserv <- regland_sp[ortho_regland[[paste0("ID.", ref_sp)]], "ratio_cons_int"]
    
    ## Correlation 
    R=cor(log2(ortho_regland[[paste0(ref_sp, "_", cell)]]+1), ortho_regland$sp_contact_conserv, method="pearson")
    rho=cor(log2(ortho_regland[[paste0(ref_sp, "_", cell)]]+1), ortho_regland$sp_contact_conserv, method="spearman")
    
    smoothScatter(log2(ortho_regland[[paste0(ref_sp, "_", cell)]]+1), ortho_regland$sp_contact_conserv, main=paste(enh, "in", cell),
                  xlab=paste(ref_sp, "gene expression log2(TPM)", sep=" "), ylab="Enhancers conserved in contact (%)", cex.lab=1.2)
    
    mtext(paste("Pearson's R = ", round(R, digits=2), ", Spearman's rho = ", round(rho, digits=2),sep=""), side=3, line=0.5, cex=0.8)
    x=log2(ortho_regland[[paste0(ref_sp, "_", cell)]]+1)
    y=ortho_regland$sp_contact_conserv
    abline(lm(y~x), col="red")
  }
}
dev.off()