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
enh = "FANTOM5"
dataset.colors=c("firebrick1", "forestgreen", "navy")
names(dataset.colors) = cells

################################################################################################################
############################################## Datas ###########################################################

ortho <- read.table(paste(pathFinalData, "SupplementaryDataset3/human2mouse_ortholog_one2one_genes_Ensembl94.txt", sep="/"), h=T, sep="\t")
rownames(ortho) <- ortho$GenestableID
ortho$dNdS <- ortho$dN/ortho$dS
ortho <- ortho[which(!is.na(ortho$dNdS) & ortho$dNdS < 50),]

expdiv=read.table(paste(path_exp, "expression_divergence/ExpressionDivergence_CellTypes_MeanTPM.txt",sep=""), h=T, stringsAsFactors=F, sep="\t")
colnames(expdiv) = c("IDHuman", "IDMouse", "adipo", "Bcell", "ESC")
expdiv <- expdiv[complete.cases(expdiv), ] # to remove if we find a solution to NA

exp_human = read.table(paste(path_exp, "human/TPMValues_CommonCellTypes.txt", sep=""), h=T, stringsAsFactors=F, row.names = 1, sep="\t")
exp_mouse = read.table(paste(path_exp, "mouse/TPMValues_CommonCellTypes.txt", sep=""), h=T, stringsAsFactors=F, row.names = 1, sep="\t")

if(ref_sp == "human"){rownames(expdiv)=expdiv$IDHuman}else{rownames(expdiv)=expdiv$IDMouse}

for (cell in cells){
  samples = grep(cell, colnames(exp_human), value=T)
  expdiv[[paste0("human_", cell, "_Mean")]] <- apply(exp_human[expdiv[,"IDHuman"], samples], 1, function(x) mean(x))
  
  if (cell == "Bcell"){cell_name="B"}else{cell_name = cell} # Because Bcell in mouse are intitulate B1, B2...
  samples = grep(cell_name, colnames(exp_mouse), value=T)
  expdiv[[paste0("mouse_", cell, "_Mean")]] <- apply(exp_mouse[expdiv[,"IDMouse"], samples], 1, function(x) mean(x))
  
  expdiv[[paste0(cell, "_Mean")]] <- (expdiv[[paste0("human_", cell, "_Mean")]] + expdiv[[paste0("mouse_", cell, "_Mean")]])/2
  
  lm1 = lm(expdiv[,cell]~log2(expdiv[[paste0(cell, "_Mean")]]+1))
  expdiv[[paste0(cell, "_ResidualExpressionDivergence")]] = lm1$residuals
  
}

################################################################################################################################
############################ PART 1 : Parallel trends among cell types  #######################################################
gene_dnds <- matrix(ncol=3, nrow=3, dimnames=list(cells, c("Mean","Conf_low", "Conf_high")))
enh_evol <- matrix(ncol=3, nrow=3, dimnames=list(cells, c("Mean","Conf_low", "Conf_high")))
contact_conserv <- matrix(ncol=3, nrow=3, dimnames=list(cells, c("Mean","Conf_low", "Conf_high")))
correl_expression <- matrix(ncol=2, nrow=3, dimnames=list(cells, c("Pearson","Spearman")))
correl_complexity <- matrix(ncol=2, nrow=3, dimnames=list(cells, c("Pearson","Spearman")))

enhancers_alignment = read.table(paste(path_evol, ref_sp, "sequence_conservation/enhancers", enh, "Alignments_stats_all_species_nonexonic_ungapped.txt", sep="/"), h=T)
enhancers_contact = read.table(paste(pathFinalData, "SupplementaryDataset4", ref_sp, enh,  "statistics_contacted_enhancers_original.txt", sep="/"), h=T, stringsAsFactors = F)
enhancers_contact$enh <-  do.call(paste,c(enhancers_contact[c("chr","start","end")],sep=":"))

for (cell in cells){
  if (cell == "Bcell"){cell_name="Bcell"}
  if (cell == "adipo"){cell_name="pre_adipo"}
  if (cell == "ESC"){cell_name="hESC"}
  
  # Enhancers contacted by gene in cell types
  enhancers_contact_in_cell <- enhancers_contact[which(enhancers_contact[[paste0(cell_name)]] > 0),]
  enh_alignment <- enhancers_alignment[which(enhancers_alignment$enh %in% enhancers_contact_in_cell$enh),]$mouse
  enh_evol[cell,] <- c(mean(enh_alignment), t.test(enh_alignment)[["conf.int"]][1], t.test(enh_alignment)[["conf.int"]][2])
  
  # dN/dS of more expressed genes 
  thirdquantile <- summary(log2(expdiv[[paste(ref_sp, cell, "Mean", sep="_")]]))[5]
  expdiv_top <- expdiv[which(log2(expdiv[[paste(ref_sp, cell, "Mean", sep="_")]]) > thirdquantile ),]
  dNdS <- ortho[which(ortho$GenestableID %in% expdiv_top$IDHuman),]$dNdS
  gene_dnds[cell,] <- c(mean(dNdS), t.test(dNdS)[["conf.int"]][1],  t.test(dNdS)[["conf.int"]][2])
  
  ################## Regulatory landscape conservation ############
  regland = read.table(paste(path_evol, ref_sp, "/evolution_summary_by_gene/", enh, "/original_evolution_summary_", cell, ".txt",sep=""), h=T, stringsAsFactors=F, sep="\t", row.names = 1)
  common=intersect(rownames(expdiv), rownames(regland))
  expdiv_cell=expdiv[common,]
  regland=regland[common,]
  
  ###### Correlation of expression ###### 
  correl_expression[cell,] = c(cor(log2(expdiv_cell[[paste(ref_sp, cell, "Mean", sep="_")]]+1), log2(expdiv_cell[[paste(target_sp, cell, "Mean", sep="_")]]+1), method="pearson"),
                               cor(log2(expdiv_cell[[paste(ref_sp, cell, "Mean", sep="_")]]+1), log2(expdiv_cell[[paste(target_sp, cell, "Mean", sep="_")]]+1), method="spearman"))
  
  ###### Calcul complexity zscore ###### 
  regland$zscore <- (regland$nb_total-mean(regland$nb_total)) / sd(regland$nb_total)
  
  # Get complexity zscore of orthologous gene in same cell
  regland_target = read.table(paste(path_evol, target_sp, "/evolution_summary_by_gene/", enh, "/original_evolution_summary_", cell, ".txt",sep=""), h=T, stringsAsFactors=F, sep="\t", row.names = 1)
  
  regland_target=regland_target[expdiv_cell$IDMouse,]
  regland_target[is.na(regland_target)] <- 0
  regland$zscore_target <- (regland_target$nb_total-mean(regland_target$nb_total)) / sd(regland_target$nb_total)
  
  ## Correlation between regulatory landscapes complexity
  correl_complexity[cell,] = c(cor(regland$zscore, regland$zscore_target, method="pearson"),
                               cor(regland$zscore, regland$zscore_target, method="spearman"))
  
  
  ######  Ratio conserved contact ###### 
  contact <- regland$nb_contact_conserv/regland$nb_total
  contact_conserv[cell,] <- c(mean(contact), t.test(contact)[["conf.int"]][1], t.test(contact)[["conf.int"]][2])
  
}

### Plot PART 1 ###
pdf(paste(pathFigures, "Figure7_bis.pdf", sep=""), width=7, height=5)
par(mai = c(0.5, 0.1, 0.5, 0.1)) # bottom, left, top, right

a <- matrix(c(1,2,3,4,5,6), nrow = 1, byrow=F)
b <- matrix(c(7,7,8,8,9,9), nrow=1, byrow=F)
c <- rbind(a,b)
layout(c)

plot.new()
legend("center", inset=c(-0.55,0), col=dataset.colors[cells], legend = cells, bty='n', lty=1, y.intersp=4)

# A - Correlation of expression level
dotchart(correl_expression[rev(cells),"Pearson"], col=dataset.colors[rev(cells)], labels='', pch=16, pt.cex=0.5,
         xlim=c(min(correl_expression)-0.02, max(correl_expression)+0.02))
mtext("Expression level \n correlation (rho)", side=1, line=3.5, cex=0.7)
mtext("A", side=3, line=1, at=0.7, font=2, cex=1)

# B - Correlattion of complexity zscore 
dotchart(correl_complexity[rev(cells),"Pearson"], col=dataset.colors[rev(cells)], labels='', pch=16, pt.cex=0.5, 
         xlim=c(min(correl_complexity)-0.02, max(correl_complexity)+0.02))
mtext("Complexity \n correlation (rho)", side=1, line=3.5, cex=0.7)
mtext("B", side=3, line=1, font=2, cex=1)

# C - dN / dS
dotchart(gene_dnds[rev(cells),"Mean"], col=dataset.colors[rev(cells)], labels='', pch=16, pt.cex=0.5, 
         xlim=c(min(gene_dnds[,"Conf_low"])-0.02, max(gene_dnds[,"Conf_high"])+0.02))
segments(x0=gene_dnds[rev(cells),"Conf_low"], x1=gene_dnds[rev(cells),"Conf_high"], y0=1:3, y1=1:3, col=dataset.colors[rev(cells)])
segments(x0=gene_dnds[rev(cells),"Conf_low"], x1=gene_dnds[rev(cells),"Conf_low"], y0=(1:3)-0.02, y1=(1:3)+0.02, col=dataset.colors[rev(cells)], lwd=1)
segments(x0=gene_dnds[rev(cells),"Conf_high"], x1=gene_dnds[rev(cells),"Conf_high"], y0=(1:3)-0.02, y1=(1:3)+0.02, col=dataset.colors[rev(cells)], lwd=1)

mtext("Top expressed gene \n dN/dS", side=1, line=3.5, cex=0.7)
mtext("C", side=3, line=1, font=2, cex=1)

# D - Enhancer Alignment
dotchart(enh_evol[rev(cells),"Mean"], col=dataset.colors[rev(cells)], pch=16, labels='', pt.cex=0.5,
         xlim=c(min(enh_evol[,"Conf_low"])-0.02, max(enh_evol[,"Conf_high"])+0.02))
segments(x0=enh_evol[rev(cells),"Conf_low"], x1=enh_evol[rev(cells),"Conf_high"], y0=1:3, y1=1:3, col=dataset.colors[rev(cells)])
segments(x0=enh_evol[rev(cells),"Conf_low"], x1=enh_evol[rev(cells),"Conf_low"], y0=(1:3)-0.02, y1=(1:3)+0.02, col=dataset.colors[rev(cells)], lwd=1)
segments(x0=enh_evol[rev(cells),"Conf_high"], x1=enh_evol[rev(cells),"Conf_high"], y0=(1:3)-0.02, y1=(1:3)+0.02, col=dataset.colors[rev(cells)], lwd=1)

mtext("Enhancer \n Alignment Score", side=1, line=3.5, cex=0.7)
mtext("D", side=3, line=1, font=2, cex=1)

# E - Conserved contacts
dotchart(contact_conserv[rev(cells),"Mean"], col=dataset.colors[rev(cells)], pch=16, labels='', pt.cex=0.5,
         xlim=c(min(contact_conserv[,"Conf_low"])-0.04, max(contact_conserv[,"Conf_high"])+0.04))
segments(x0=contact_conserv[rev(cells),"Conf_low"], x1=contact_conserv[rev(cells),"Conf_high"], y0=1:3, y1=1:3, col=dataset.colors[rev(cells)], lwd=1)
segments(x0=contact_conserv[rev(cells),"Conf_low"], x1=contact_conserv[rev(cells),"Conf_low"], y0=(1:3)-0.02, y1=(1:3)+0.02, col=dataset.colors[rev(cells)], lwd=1)
segments(x0=contact_conserv[rev(cells),"Conf_high"], x1=contact_conserv[rev(cells),"Conf_high"], y0=(1:3)-0.02, y1=(1:3)+0.02, col=dataset.colors[rev(cells)], lwd=1)

mtext("Proportion of \n conserved contact", side=1, line=3.5, cex=0.7)
mtext("E", side=3, line=1, font=2, cex=1)


################################################################################################################################
############################ PART2 : Expression divergence vs Regulatory Landscape divergence  #################################
################################################################################################################################
par(mai = c(0.5, 0.5, 0.5, 0.1)) # bottom, left, top, right

data_cell <- list()

for (cell in cells){
  regland = read.table(paste(path_evol, ref_sp, "/evolution_summary_by_gene/", enh, "/original_evolution_summary_", cell, ".txt",sep=""), h=T, stringsAsFactors=F, sep="\t", row.names = 1)
  common=intersect(rownames(expdiv), rownames(regland))
  expdiv_cell=expdiv[common,]
  regland=regland[common,]
  
  # Made decile of nb enhancers
  regland$class_nb_contact <- cut2(regland$nb_total, g=5)

  regland$ratio_cons_seq = regland$nb_seq_conserv/regland$nb_total
  regland$ratio_cons_int = ifelse(regland$nb_seq_conserv > 0, regland$nb_contact_conserv/regland$nb_seq_conserv, 0)
  
  regland$class_cons_seq=cut(regland$ratio_cons_seq, breaks=c(0, 0.001, 0.25, 0.50, 0.75, 1), include.lowest=T)
  regland$class_cons_int=cut(regland$ratio_cons_int,  breaks=c(0, 0.001, 0.25, 0.50, 0.75, 1), include.lowest=T)
  
  regland$divergence <- expdiv_cell[[paste0(cell)]]
  
  data_cell[[cell]] <- regland
}

######################################################################################################################
#################################### F - Expression Divergence vs Number of enhancers ################################
mean_divergence <- t(sapply(data_cell, function(x)   tapply(x$divergence, as.factor(x$class_nb_contact), mean, na.rm=T)))
divergence_conf_low <- t(sapply(data_cell, function(x) tapply(x$divergence, as.factor(x$class_nb_contact), function(y) {z<-t.test(y); return(z[["conf.int"]][1])})))
divergence_conf_high <- t(sapply(data_cell, function(x) tapply(x$divergence, as.factor(x$class_nb_contact), function(y) {z<-t.test(y); return(z[["conf.int"]][2])})))

plot(as.numeric(mean_divergence["Bcell",]), type="l", col=dataset.colors["Bcell"], ylim=c(0.5,0.75), xlab="", ylab="", axes=F)
lines(as.numeric(mean_divergence["ESC",]), col=dataset.colors["ESC"], lwd=1.5)
lines(as.numeric(mean_divergence["adipo",]), col=dataset.colors["adipo"], lwd=1.5)

## X axis
axis(side=1, at=1:5, labels=1:5, mgp=c(3, 0.65, 0), cex.axis=0.8)
axis(side=2, mgp=c(3, 0.75, 0), cex.axis=0.8)

## axis labels
mtext("Quantile of Complexity", side=1, line=2.5, cex=0.7)
mtext("Expression Divergence", side=2, line=2.7, cex=0.7)

## confidence intervals

for(dataset in rownames(mean_divergence)){
  for (col in 1:ncol(mean_divergence)){
    segments(x0=col, x1=col, y0=as.numeric(divergence_conf_low[dataset,col]), y1=as.numeric(divergence_conf_high[dataset,col]), col=dataset.colors[dataset])
  }
}

## legend & plot label
#legend("topright", legend=names(dataset.colors), col=dataset.colors, lty=1, bty='n', inset=c(0.05, 0), xpd=NA, cex=1.4)
mtext("F", side=3, line=1, at=1, font=2, cex=1.2)

######################################################################################################################
################################  G - Expression Divergence vs Number of conserved enhancers ######################## 
mean_divergence <- t(sapply(data_cell, function(x)   tapply(x$divergence, as.factor(x$class_cons_seq), mean, na.rm=T)))
divergence_conf_low <- t(sapply(data_cell, function(x) tapply(x$divergence, as.factor(x$class_cons_seq), function(y) {z<-t.test(y); return(z[["conf.int"]][1])})))
divergence_conf_high <- t(sapply(data_cell, function(x) tapply(x$divergence, as.factor(x$class_cons_seq), function(y) {z<-t.test(y); return(z[["conf.int"]][2])})))

plot(as.numeric(mean_divergence["Bcell",]), type="l", col=dataset.colors["Bcell"], ylim=c(0.5,0.75), xlab="", ylab="", axes=F)
lines(as.numeric(mean_divergence["ESC",]), col=dataset.colors["ESC"], lwd=1.5)
lines(as.numeric(mean_divergence["adipo",]), col=dataset.colors["adipo"], lwd=1.5)

## X axis
breaks_names=c(">0.1", "0.1-25", "26-50", "51-75", ">75")
axis(side=1, at=1:5, labels=breaks_names, mgp=c(3, 0.65, 0), cex.axis=0.8)
axis(side=2, mgp=c(3, 0.75, 0), cex.axis=0.8)

## axis labels
mtext("Sequences conserved (%)", side=1, line=2.25, cex=0.7)

## confidence intervals

for(dataset in rownames(mean_divergence)){
  for (col in 1:ncol(mean_divergence)){
    segments(x0=col, x1=col, y0=as.numeric(divergence_conf_low[dataset,col]), y1=as.numeric(divergence_conf_high[dataset,col]), col=dataset.colors[dataset])
  }
}

## plot label
mtext("G", side=3, line=1, at=1, font=2, cex=1)

######################################################################################################################
############################## H - Expression Divergence vs Number of conserved contacts ############################# 
mean_divergence <- t(sapply(data_cell, function(x)   tapply(x$divergence, as.factor(x$class_cons_int), mean, na.rm=T)))
divergence_conf_low <- t(sapply(data_cell, function(x) tapply(x$divergence, as.factor(x$class_cons_int), function(y) {z<-t.test(y); return(z[["conf.int"]][1])})))
divergence_conf_high <- t(sapply(data_cell, function(x) tapply(x$divergence, as.factor(x$class_cons_int), function(y) {z<-t.test(y); return(z[["conf.int"]][2])})))

plot(as.numeric(mean_divergence["Bcell",]), type="l", col=dataset.colors["Bcell"], ylim=c(0.45,0.7), xlab="", ylab="", axes=F)
lines(as.numeric(mean_divergence["ESC",]), col=dataset.colors["ESC"], lwd=1.5)
lines(as.numeric(mean_divergence["adipo",]), col=dataset.colors["adipo"], lwd=1.5)

## X axis
breaks_names=c(">0.1", "0.1-25", "26-50", "51-75", ">75")
axis(side=1, at=1:5, labels=breaks_names, mgp=c(3, 0.65, 0), cex.axis=0.8)
axis(side=2, mgp=c(3, 0.75, 0), cex.axis=0.8)

## axis labels
mtext("Contacts conserved (%)", side=1, line=2.25, cex=0.7)

## confidence intervals

for(dataset in rownames(mean_divergence)){
  for (col in 1:ncol(mean_divergence)){
    segments(x0=col, x1=col, y0=as.numeric(divergence_conf_low[dataset,col]), y1=as.numeric(divergence_conf_high[dataset,col]), col=dataset.colors[dataset])
  }
}

## plot label
mtext("H", side=3, line=1, at=1, font=2, cex=1)

dev.off()
