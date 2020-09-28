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

ortho <- read.table(paste(pathFinalData, "SupplementaryDataset7/human/gene_orthology/human2mouse_orthologue_dNdS.txt", sep="/"), h=T, sep="\t")
rownames(ortho) <- ortho$GenestableID
ortho$dNdS <- ortho$dN/ortho$dS
ortho <- ortho[which(!is.na(ortho$dNdS) & ortho$dNdS < 50),]

expdiv=read.table(paste(path_exp, "expression_divergence/ExpressionDivergence_CellTypes_MedianTPM.txt",sep=""), h=T, stringsAsFactors=F, sep="\t")
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
  if (cell == "Bcell"){if (ref_sp == "human"){cell_name=c("Bcell", "TB", "NB")}else{cell_name=c("preB_aged", "preB_young")}}
  if (cell == "adipo"){if (ref_sp == "human"){cell_name=c("pre_adipo")}else{cell_name=c("preadip_4H", "preadip_D0", "preadip_D2")}}
  if (cell == "ESC"){if (ref_sp == "human"){cell_name=c("hESC")}else{cell_name=c("ESC", "ESC_18", "ESC_NKO")}}
  
  # Enhancers contacted by gene in cell types
  enhancers_contact_in_cell <- enhancers_contact[which(apply(as.matrix(enhancers_contact[,cell_name]), 1, function(X) any(X > 0))),]
  enh_alignment <- enhancers_alignment[which(enhancers_alignment$enh %in% enhancers_contact_in_cell$enh),target_sp]
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
  
  if (target_sp == "mouse"){ID = "IDMouse"}else{ID = "IDHuman"}
  
  regland_target=regland_target[expdiv_cell[[ID]],]
  regland_target[is.na(regland_target)] <- 0
  regland$zscore_target <- (regland_target$nb_total-mean(regland_target$nb_total)) / sd(regland_target$nb_total)
  
  ## Correlation between regulatory landscapes complexity
  correl_complexity[cell,] = c(cor(regland$zscore, regland$zscore_target, method="pearson"),
                               cor(regland$zscore, regland$zscore_target, method="spearman"))
  
  ######  Ratio conserved contact ###### 
  regland <- regland[which(regland$nb_total > 5),]
  contact <- regland$nb_contact_conserv/regland$nb_total
  contact_conserv[cell,] <- c(mean(contact), t.test(contact)[["conf.int"]][1], t.test(contact)[["conf.int"]][2])
  
}

### Plot PART 1 ###
#pdf(paste(pathFigures, "Figure7.pdf", sep=""), width=7, height=5)
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

# B - Correlation of complexity zscore 
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
par(mfrow=c(1,3))
for (cell in cells){
  #file = paste(path_evol, ref_sp, "/evolution_summary_by_gene/", enh, "/original_evolution_summary_", cell, ".txt", sep="")
  file = paste("/home/laverre/Data/Regulatory_landscape/result/Supplementary_dataset6_regulatory_landscape_evolution/human/", enh, "_original_evolution_summary_", cell, "_100kb.txt", sep="")
  regland = read.table(file, h=T, stringsAsFactors=F, sep="\t", row.names = 1)

  #firstquantile <- summary(log2(expdiv[[paste(ref_sp, cell, "Mean", sep="_")]]))[1]
  #expdiv_cell <- expdiv[which(log2(expdiv[[paste(ref_sp, cell, "Mean", sep="_")]]) > firstquantile ),]
  
  common=intersect(rownames(expdiv_cell), rownames(regland))
  expdiv_cell=expdiv_cell[common,]
  regland=regland[common,]
  

  # Made decile of nb enhancers
  if (enh == "ENCODE"){BREAKS=c(1, 5, 10, 20, 30, max(regland$nb_total))}else{BREAKS=c(1, 5, 8, 10, 15, max(regland$nb_total))}
  
  regland$class_nb_contact <- cut(regland$nb_total, breaks=BREAKS, include.lowest=T)
  regland$class_distance <- cut(regland$median_dist, breaks=c(0, 100000, 250000, 500000, 1000000, max(regland$median_dist)), include.lowest=T)
  regland$divergence <- expdiv_cell[[paste0(cell)]] #, "_ResidualExpressionDivergence"
  regland$class_expression <- cut2(expdiv_cell[[paste0("human_", cell, "_Mean")]], g=3)

  #regland <- regland[which(regland$nb_total > 5),]

  regland$ratio_cons_seq = regland$nb_seq_conserv/regland$nb_total
  regland$ratio_cons_synt = ifelse(regland$nb_seq_conserv > 0, regland$nb_synt2M_conserv/regland$nb_seq_conserv, 0)
  regland$ratio_cons_int = ifelse(regland$nb_seq_conserv > 0, regland$nb_contact_conserv/regland$nb_seq_conserv, 0)
  
  regland$class_cons_seq=cut(regland$ratio_cons_seq, breaks=c(0, 0.25, 0.4, 0.50, 0.75, 1), include.lowest=T)
  regland$class_cons_synt=cut(regland$ratio_cons_synt, breaks=c(0, 0.75, 1), include.lowest=T)
  regland$class_cons_int=cut(regland$ratio_cons_int,  breaks=c(0, 0.01, 0.25, 0.5, 0.75, 1), include.lowest=T)
  
  data_cell[[cell]] <- regland
  
}


######################################################################################################################
#################################### F - Expression Divergence vs Number of enhancers ################################
mean_divergence <- t(sapply(data_cell, function(x)   tapply(x$divergence, as.factor(x$class_nb_contact), mean, na.rm=T)))
divergence_conf_low <- t(sapply(data_cell, function(x) tapply(x$divergence, as.factor(x$class_nb_contact), function(y) {z<-t.test(y); return(z[["conf.int"]][1])})))
divergence_conf_high <- t(sapply(data_cell, function(x) tapply(x$divergence, as.factor(x$class_nb_contact), function(y) {z<-t.test(y); return(z[["conf.int"]][2])})))

if (enh == "ENCODE"){YLIM=c(0.4,0.7)}else{YLIM=c(0.2,0.75)}
#if (enh == "ENCODE"){YLIM=c(-0.05,0.05)}else{YLIM=c(0.2,0.75)}
plot(as.numeric(mean_divergence["Bcell",]), type="l", col=dataset.colors["Bcell"], ylim=YLIM, xlab="", ylab="", axes=F)
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

plot(as.numeric(mean_divergence["Bcell",]), type="l", col=dataset.colors["Bcell"], ylim=YLIM, xlab="", ylab="", axes=F)
lines(as.numeric(mean_divergence["ESC",]), col=dataset.colors["ESC"], lwd=1.5)
lines(as.numeric(mean_divergence["adipo",]), col=dataset.colors["adipo"], lwd=1.5)

## X axis
breaks_names=c(">25", "25-40", "41-50", "51-75", ">75")

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

plot(as.numeric(mean_divergence["Bcell",]), type="l", col=dataset.colors["Bcell"], ylim=YLIM, xlab="", ylab="", axes=F)
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

#dev.off()

######################################################################################################################
########################################## Cardoso - Moreira Divergence ##########################################
setwd("/home/laverre/Manuscript/SupplementaryDataset6/expression_divergence")

expdiv=read.table("ExpressionDivergence_CardosoMoreira2019_correlations.txt", h=T, stringsAsFactors=F, sep="\t")
rownames(expdiv)=expdiv$IDHuman
enh="ENCODE"
maxdist = "25000"
file = paste("/home/laverre/Data/Regulatory_landscape/result/Supplementary_dataset6_regulatory_landscape_evolution/human/", enh, "_original_evolution_summary_all_", maxdist, ".txt", sep="")
regland = read.table(file, h=T, stringsAsFactors=F, sep="\t", row.names = 1)

common=intersect(rownames(expdiv), rownames(regland))
expdiv=expdiv[common,]
regland=regland[common,]

regland$ratio_cons_seq = regland$nb_seq_conserv/regland$nb_total
regland$ratio_cons_synt = ifelse(regland$nb_seq_conserv > 0, regland$nb_synt2M_conserv/regland$nb_seq_conserv, NA)
regland$ratio_cons_int = ifelse(regland$nb_seq_conserv > 0, regland$nb_contact_conserv/regland$nb_seq_conserv, 0)

regland$class_nb_contact=cut(regland$nb_total, breaks=c(1,5,10,25,50,75, max(regland$nb_total)), include.lowest=T)
regland$class_cons_seq=cut(regland$ratio_cons_seq, breaks=c(0, 0.25, 0.4, 0.50, 0.75, 1), include.lowest=T)
regland$class_cons_synt=cut(regland$ratio_cons_synt, breaks=c(0, 0.75, 1), include.lowest=T)
regland$class_cons_int=cut(regland$ratio_cons_int,  breaks=c(0, 0.01, 0.25, 0.5, 0.75, 1), include.lowest=T)

# Divergence ~ Nb contact
par(mfrow=c(1,1))
breaks_names = c("1-5","6-10","11-25", "26-50", "51-75",">75")

boxplot(expdiv$CorrelationSpearman~regland$class_nb_contact, outline=F, notch=T,
        xlab="Total number of contacts", ylab="Spearman's correlation", names=breaks_names, boxwex = 0.7)

boxplot(expdiv$ResidualSpearman~regland$class_nb_contact, outline=F, notch=T, boxwex = 0.7,
        xlab="Total number of contacts", ylab="Residual expression divergence", names=breaks_names)

# Divergence ~ Nb seq conserv
regland <- regland[which(regland$nb_total >= 5),]
common=intersect(rownames(expdiv), rownames(regland))
expdiv=expdiv[common,]
regland=regland[common,]

par(mfrow=c(1,2))
boxplot(expdiv$CorrelationSpearman~regland$class_cons_seq, outline=F, notch=T,
        xlab="Ratio of conserved sequences", ylab="Spearman's correlation", boxwex = 0.7)

boxplot(expdiv$ResidualSpearman~regland$class_cons_seq, outline=F, notch=T, boxwex = 0.7,
        xlab="Ratio of conserved sequences", ylab="Residual Spearman's correlation")

# Divergence ~ Nb synteny conserv
par(mfrow=c(1,2))
boxplot(expdiv$CorrelationSpearman~regland$class_cons_synt, outline=F, notch=T,
        xlab="Ratio of conserved synteny", ylab="Spearman's correlation", boxwex = 0.7)

boxplot(expdiv$ResidualExpressionDivergence~regland$class_cons_synt, outline=F, notch=T, boxwex = 0.7,
        xlab="Ratio of conserved synteny", ylab="Residual Spearman's correlation")

# Divergence ~ Nb contact conserv
regland <- regland[which(regland$nb_seq_conserv >= 5),]
common=intersect(rownames(expdiv), rownames(regland))
expdiv=expdiv[common,]
regland=regland[common,]

par(mfrow=c(1,2))

boxplot(expdiv$CorrelationSpearman~regland$class_cons_int, outline=F, notch=T,
        xlab="Ratio of conserved contact", ylab="Spearman's correlation", boxwex = 0.7)

boxplot(expdiv$ResidualSpearman~regland$class_cons_int, outline=F, notch=T, boxwex = 0.7,
        xlab="Ratio of conserved contact", ylab="Residual Spearman's correlation")

write.table(regland[order(-regland$ratio_cons_int),11:13],"test.txt",sep="\t",row.names=TRUE, quote=F)
write.table(expdiv[order(-expdiv$ExpressionDivergence),3:4],"test_div_expression.txt",sep="\t",row.names=TRUE, quote=F)
write.table(expdiv[order(-expdiv$Bcell),3:4],"test_div_expression_Bcell.txt",sep="\t",row.names=TRUE, quote=F)


plot(expdiv$ExpressionDivergence, expdiv$CorrelationSpearman, cex=0.6)
rho=cor(expdiv$ExpressionDivergence, expdiv$CorrelationSpearman, method="spearman")
R=cor(expdiv$ExpressionDivergence, expdiv$CorrelationSpearman, method="pearson")

mtext(paste("R2 = ", round(R, digits=2), ", rho = ",round(rho, digits=2),sep=""), side=3, line=0.5)
abline(lm(expdiv$ExpressionDivergence~expdiv$CorrelationSpearman), col="red")
