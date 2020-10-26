######################################################################################################################
library(Hmisc)

#setwd("/home/laverre/Data/Regulatory_landscape/scripts/main_figures")
source("parameters.R") ## pathFinalData are defined based on the user name

sp="human"

load(paste(pathFigures, "RData/data.gene.annotations.RData", sep=""))
load(paste(pathFigures, "RData/Fig6_", sp, "_common_cells.Rdata", sep=""))
load(paste(pathFigures, "RData/Fig6_", sp, "_all_samples_CardosoMoreira.Rdata", sep=""))

if (sp == "human"){sp_name="Human"}else{sp_name="Mouse"}
DivergenceMeasure = "CorrelationSpearman" # "EuclidianSimilarity" or "CorrelationSpearman" or SpearmanResidual

cells <- c("Bcell", "ESC", "adipo")
dataset.colors=c("firebrick1", "forestgreen", "navy")
names(dataset.colors) = cells

################################################################################################################################
#pdf(file=paste(pathFigures, "Figure7.pdf", sep=""))

m=matrix(rep(NA, 8*12), nrow=8)
m[1,]=c(rep(1, 5), rep(3,3), rep(4,4))
m[2,]=c(rep(2, 5), rep(3,3), rep(4,4))
m[3,]=c(rep(5, 3), rep(6,3),  rep(7,3),  rep(8,3))
m[4,]=c(rep(5, 3), rep(6,3),  rep(7,3),  rep(8,3))
m[5,]=c(rep(9, 2), rep(10,2),  rep(11,2),  rep(12,2), rep(13,2), rep(14,2))
m[6,]=c(rep(9, 2), rep(10,2),  rep(11,2),  rep(12,2), rep(13,2), rep(14,2))
m[7,]=c(rep(15, 4), rep(16,4),  rep(17,4))
m[8,]=c(rep(15, 4), rep(16,4),  rep(17,4))
layout(m)

################################################################################################################################
############################## PART1 : All cells & Cardoso-Moreira  ##########################################################

## A - Relative expression profile SHH 
tissues.human=unlist(lapply(samples.human, function(x) unlist(strsplit(x, split="_"))[2]))
tissues.mouse=unlist(lapply(samples.mouse, function(x) unlist(strsplit(x, split="_"))[2]))
stage.human=unlist(lapply(samples.human, function(x) unlist(strsplit(x, split="_"))[3]))
stage.mouse=unlist(lapply(samples.mouse, function(x) unlist(strsplit(x, split="_"))[3]))

col.tissues=c("navy", "steelblue", "indianred", "seagreen", "orange")
names(col.tissues)=c("Brain", "Cerebellum", "Heart", "Kidney", "Liver")

rel.shh.human=shh.human/sum(shh.human)
rel.shh.mouse=shh.mouse/sum(shh.mouse)

par(mar=c(2,4.5,2.5,1.1)) # bottom, left, top, right
b=barplot(rel.shh.human, col=col.tissues[tissues.human], border=col.tissues[tissues.human], axes=F)
axis(side=2, cex=0.95)
mtext("Human", side=2, line=3, cex=1)
mtext(stage.human[seq(1,length(stage.human),2)], side=1, at=b[seq(1,length(b),2)], las=2, line=0.5, cex=0.6)
mtext("A", side=3, at=0.5, font=2, cex=1.1, line=1.1)
legend("topright", legend=c("forebrain", "cerebellum", "kidney", "liver", "heart"), ncol=2,
       fill=col.tissues[c("Brain", "Cerebellum",  "Kidney", "Liver","Heart")],
       border=col.tissues[c("Brain", "Cerebellum",  "Kidney", "Liver","Heart")], inset=c(-0.03, -0.11), bty="n")
mtext("SHH expression", side=3, line=1, cex=0.8)

b=barplot(rel.shh.mouse, col=col.tissues[tissues.mouse], border=col.tissues[tissues.mouse], axes=F)
axis(side=2, cex=0.95)
mtext("Mouse ", side=2, line=3, cex=1)
mtext(stage.mouse[seq(1,length(stage.mouse),2)], side=1, at=b[seq(1,length(b),2)], las=2, line=0.5, cex=0.6)

## B - classes of expression specificity 
xpos=seq(1, length(levels(expdiv_all$classTau)), 1)
names(xpos) = levels(expdiv_all$classTau)

smallx=c(-0.15, -0.075, 0.075, 0.15)
names(smallx)=enhancer.datasets[[sp]]

if (DivergenceMeasure == "CorrelationSpearman"){ylim=c(-0.1, 1)}else{ylim=c(0.7, 1.2)}
xlim=c(0.5, length(levels(expdiv_all$classTau))+0.5)

par(mar=c(4.5, 3.75, 3, 1.1))
plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=ylim, xaxs="i", yaxs="i")

for(class in levels(expdiv_all$classTau)){
  x=xpos[class]
  boxplot(expdiv_all[which(expdiv_all$classTau == class), DivergenceMeasure], at=x, axes=F, add=T, notch=T, outline=F)
}

abline(v=xpos[1:length(levels(expdiv_all$classTau))-1]+0.5, lty=3, col="gray40")

axis(side=1, at=xpos, mgp=c(3, 0.5, 0), labels=rep("", length(levels(expdiv_all$classTau))), cex.axis=0.8)
mtext(c("generalist", "", "", "specialist"), at=xpos, side=1, line=1, cex=0.7)
mtext("Expression Specificity", side=1, line=2.5, cex=0.9)

axis(side=2, mgp=c(3, 0.75, 0), cex.axis=1.1)
mtext("Expression Similarity", side=2, line=2.5, cex=0.9)

mtext("B", side=3, at=0.5, font=2, cex=1.1, line=1.5)

## C - Specificity vs nb enhancers
ylim=c(0, 35)
par(mar=c(4.5, 3.75, 3, 1.1))
plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=ylim, xaxs="i", yaxs="i")

for(enh in enhancer.datasets[[sp]]){
  
  for(class in levels(expdiv[[enh]]$classTau)){
    this.genes=rownames(expdiv[[enh]][which(expdiv[[enh]]$classTau == class),])
    
    x=xpos[class]+smallx[enh]
    
    b=boxplot(regland[[enh]][this.genes, "nb_total"], plot=FALSE)
    med=median(regland[[enh]][this.genes, "nb_total"])
    ci=as.numeric(b$conf)
    
    points(x, med, pch=20, col=col.enhancers[enh], cex=1.1)
    segments(x, ci[1], x, ci[2], col=col.enhancers[enh])
  }
}

abline(v=xpos[1:length(levels(expdiv_all$classTau))-1]+0.5, lty=3, col="gray40")
axis(side=1, at=xpos, mgp=c(3, 0.5, 0), labels=rep("", length(levels(expdiv_all$classTau))), cex.axis=0.8)
mtext(c("generalist", "", "", "specialist"), at=xpos, side=1, line=1, cex=0.7)
mtext("Expression Specificity", side=1, line=2.5, cex=0.9)

axis(side=2, mgp=c(3, 0.75, 0), cex.axis=1.1)
mtext("nb. contacted enhancers", side=2, line=2.5, cex=0.9)

legend("topright", legend=label.enhancers[enhancer.datasets[[sp]]], lty=1, 
       col=col.enhancers[enhancer.datasets[[sp]]], bty="o", cex=1.1, seg.len=1,
       inset=c(-0.01, -0.2), xpd=NA, box.col="white", bg="white")
mtext("C", side=3, at=0.5, font=2, cex=1.1, line=1.5)
#, inset=c(-0.01, -0.25)

## D - Similarity vs nb enhancers
par(mai = c(0.5, 0.5, 0.5, 0.1)) # bottom, left, top, right

if (DivergenceMeasure == "CorrelationSpearman"){ylim=c(0.45, 0.65)}else{ylim=c(-0.03, 0.1)}

xlim=c(0.5, length(levels(regland[[enh]]$class_nb_contact))+0.5)

plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=ylim, xaxs="i", yaxs="i")

for(enh in enhancer.datasets[[sp]]){
  message(enh)
  for(class in levels(regland[[enh]]$class_nb_contact)){
    this.genes=rownames(regland[[enh]][which(regland[[enh]]$class_nb_contact == class),])
    
    xpos=seq(1, length(levels(regland[[enh]]$class_nb_contact)), 1)
    names(xpos) = levels(regland[[enh]]$class_nb_contact)
    x=xpos[class]+smallx[enh]
    
    b=boxplot(expdiv[[enh]][this.genes, DivergenceMeasure], plot=FALSE)
    med=median(expdiv[[enh]][this.genes, DivergenceMeasure])
    ci=as.numeric(b$conf)
    
    points(x, med, pch=20, col=col.enhancers[enh], cex=1.1)
    segments(x, ci[1], x, ci[2],  col=col.enhancers[enh])
  }
}

abline(v=xpos[1:4]+0.5, lty=3, col="gray40")
axis(side=2, mgp=c(3, 0.75, 0), cex.axis=1.1)
mtext("Spearman's coefficient", side=2, line=2.5, cex=0.9)
axis(side=1, at=xpos, mgp=c(3, 0.5, 0), labels=rep("", length(levels(regland[[enh]]$class_nb_contact))), cex.axis=0.8)
mtext(xpos, at=xpos, side=1, line=1, cex=0.8)
mtext("Complexity", side=1, line=2.5, cex=0.9)
mtext("D", side=3, at=0.45, font=2, cex=1.1, line=0.5)

## E - Similarity vs enhancers conserved in sequences 
xlim=c(0.5, length(levels(regland[[enh]]$class_align_score))+0.5)
par(mai = c(0.5, 0.1, 0.5, 0.1)) # bottom, left, top, right

plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=ylim, xaxs="i", yaxs="i")

for(enh in enhancer.datasets[[sp]]){
  message(enh)
  for(class in levels(regland[[enh]]$class_align_score)){
    this.genes=rownames(regland[[enh]][which(regland[[enh]]$class_align_score == class),])
    
    xpos=seq(1, length(levels(regland[[enh]]$class_align_score)), 1)
    names(xpos) = levels(regland[[enh]]$class_align_score)
    x=xpos[class]+smallx[enh]
    
    b=boxplot(expdiv[[enh]][this.genes, DivergenceMeasure], plot=FALSE)
    med=median(expdiv[[enh]][this.genes, DivergenceMeasure])
    ci=as.numeric(b$conf)
    
    points(x, med, pch=20, col=col.enhancers[enh], cex=1.1)
    segments(x, ci[1], x, ci[2], col=col.enhancers[enh])
  }
}

abline(v=xpos[1:4]+0.5, lty=3, col="gray40")
axis(side=1, at=xpos, mgp=c(3, 0.5, 0), labels=rep("", length(levels(regland[[enh]]$class_align_score))), cex.axis=0.8)
mtext(xpos, at=xpos, side=1, line=1, cex=0.8)
mtext("Alignment score", side=1, line=2.5, cex=0.9)
mtext("E", side=3, at=0.45, font=2, cex=1.1, line=0.5)

## F - Similarity vs enhancers conserved in synteny
xlim=c(0.5, length(levels(regland[[enh]]$class_cons_synt))+0.5)

plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=ylim, xaxs="i", yaxs="i")

for(enh in enhancer.datasets[[sp]]){
  message(enh)
  for(class in levels(regland[[enh]]$class_cons_synt)){
    this.genes=rownames(regland[[enh]][which(regland[[enh]]$class_cons_synt == class),])
    
    xpos=seq(1, length(levels(regland[[enh]]$class_cons_synt)), 1)
    names(xpos) = levels(regland[[enh]]$class_cons_synt)
    x=xpos[class]+smallx[enh]
    
    b=boxplot(expdiv[[enh]][this.genes, DivergenceMeasure], plot=FALSE)
    med=median(expdiv[[enh]][this.genes, DivergenceMeasure])
    ci=as.numeric(b$conf)
    
    points(x, med, pch=20, col=col.enhancers[enh], cex=1.1)
    segments(x, ci[1], x, ci[2], col=col.enhancers[enh])
  }
}

abline(v=xpos[1:2]+0.5, lty=3, col="gray40")
axis(side=1, at=xpos, mgp=c(3, 0.5, 0), labels=rep("", length(levels(regland[[enh]]$class_cons_synt))), cex.axis=0.8)
mtext(c("<75%", "75-99%", ">99%"), at=xpos, side=1, line=1, cex=0.8)
mtext("Synteny Conservation", side=1, line=2.5, cex=0.9)
mtext("F", side=3, at=0.45, font=2, cex=1.1, line=0.5)

### G - Similarity vs enhancers conserved in contact 
xlim=c(0.5, length(levels(regland[[enh]]$class_cons_int))+0.5)
plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=ylim, xaxs="i", yaxs="i")

for(enh in enhancer.datasets[[sp]]){
  for(class in levels(regland[[enh]]$class_cons_int)){
    this.genes=rownames(regland[[enh]][which(regland[[enh]]$class_cons_int == class),])
    
    xpos=seq(1, length(levels(regland[[enh]]$class_cons_int)), 1)
    names(xpos) = levels(regland[[enh]]$class_cons_int)
    x=xpos[class]+smallx[enh]
    
    b=boxplot(expdiv[[enh]][this.genes, DivergenceMeasure], plot=FALSE)
    med=median(expdiv[[enh]][this.genes, DivergenceMeasure])
    ci=as.numeric(b$conf)
    
    points(x, med, pch=20, col=col.enhancers[enh], cex=1.1)
    segments(x, ci[1], x, ci[2], col=col.enhancers[enh])
  }
}

abline(v=xpos[1:4]+0.5, lty=3, col="gray40")
axis(side=1, at=xpos, mgp=c(3, 0.5, 0), labels=rep("", length(levels(regland[[enh]]$class_cons_int))), cex.axis=0.8)
mtext(c("<1%", "1-25%", "25-50%", "50-75%", ">75%"), at=xpos, side=1, line=1, cex=0.8)
mtext("Contact Conservation", side=1, line=2.5, cex=0.9)
mtext("G", side=3, at=0.45, font=2, cex=1.1, line=0.5)


################################################################################################################################
######################################### PART2 : Common cells types  ##########################################################

par(mai = c(0.5, 0.1, 0.5, 0.1)) # bottom, left, top, right

#### H - Parallel trends among cell types
plot.new()
legend("center",col=dataset.colors[cells], legend = cells, bty='n', lty=1, y.intersp=4, cex=1.2)
mtext("H", side=3, line=1, at=0.7, font=2, cex=1)

#  Correlation of expression level
dotchart(correl_expression[rev(cells),"Pearson"], col=dataset.colors[rev(cells)], labels='', pch=16, pt.cex=0.8,
         xlim=c(min(correl_expression)-0.02, max(correl_expression)+0.02))
mtext("Expression level \n correlation (rho)", side=1, line=3.5, cex=0.8)

#  Correlattion of complexity zscore 
dotchart(correl_complexity[rev(cells),"Pearson"], col=dataset.colors[rev(cells)], labels='', pch=16, pt.cex=0.8, 
         xlim=c(min(correl_complexity), max(correl_complexity)+0.01))
mtext("Complexity \n correlation (rho)", side=1, line=3.5, cex=0.8)

#  dN / dS
dotchart(gene_dnds[rev(cells),"Mean"], col=dataset.colors[rev(cells)], labels='', pch=16, pt.cex=0.8, 
         xlim=c(min(gene_dnds[,"Conf_low"])-0.02, max(gene_dnds[,"Conf_high"])+0.02))
segments(x0=gene_dnds[rev(cells),"Conf_low"], x1=gene_dnds[rev(cells),"Conf_high"], y0=1:3, y1=1:3, col=dataset.colors[rev(cells)])
segments(x0=gene_dnds[rev(cells),"Conf_low"], x1=gene_dnds[rev(cells),"Conf_low"], y0=(1:3)-0.02, y1=(1:3)+0.02, col=dataset.colors[rev(cells)], lwd=1)
segments(x0=gene_dnds[rev(cells),"Conf_high"], x1=gene_dnds[rev(cells),"Conf_high"], y0=(1:3)-0.02, y1=(1:3)+0.02, col=dataset.colors[rev(cells)], lwd=1)

mtext("dN/dS", side=1, line=2.5, cex=0.8)

#  Enhancer Alignment
dotchart(enh_evol[rev(cells),"Mean"], col=dataset.colors[rev(cells)], pch=16, labels='', pt.cex=0.8,
         xlim=c(min(enh_evol[,"Conf_low"])-0.02, max(enh_evol[,"Conf_high"])+0.02))
segments(x0=enh_evol[rev(cells),"Conf_low"], x1=enh_evol[rev(cells),"Conf_high"], y0=1:3, y1=1:3, col=dataset.colors[rev(cells)])
segments(x0=enh_evol[rev(cells),"Conf_low"], x1=enh_evol[rev(cells),"Conf_low"], y0=(1:3)-0.02, y1=(1:3)+0.02, col=dataset.colors[rev(cells)], lwd=1)
segments(x0=enh_evol[rev(cells),"Conf_high"], x1=enh_evol[rev(cells),"Conf_high"], y0=(1:3)-0.02, y1=(1:3)+0.02, col=dataset.colors[rev(cells)], lwd=1)

mtext("Alignment Score", side=1, line=2.5, cex=0.8)

#  Conserved contacts
dotchart(contact_conserv[rev(cells),"Mean"], col=dataset.colors[rev(cells)], pch=16, labels='', pt.cex=0.8,
         xlim=c(min(contact_conserv[,"Conf_low"])-0.04, max(contact_conserv[,"Conf_high"])+0.04))
segments(x0=contact_conserv[rev(cells),"Conf_low"], x1=contact_conserv[rev(cells),"Conf_high"], y0=1:3, y1=1:3, col=dataset.colors[rev(cells)], lwd=1)
segments(x0=contact_conserv[rev(cells),"Conf_low"], x1=contact_conserv[rev(cells),"Conf_low"], y0=(1:3)-0.02, y1=(1:3)+0.02, col=dataset.colors[rev(cells)], lwd=1)
segments(x0=contact_conserv[rev(cells),"Conf_high"], x1=contact_conserv[rev(cells),"Conf_high"], y0=(1:3)-0.02, y1=(1:3)+0.02, col=dataset.colors[rev(cells)], lwd=1)

mtext("Conserved contact", side=1, line=2.5, cex=0.8)


par(mai = c(0.5, 0.5, 0.5, 0.1)) # bottom, left, top, right
######################## F - Expression Divergence vs Number of enhancers 
divergence="residual"
if (divergence == "residual"){YLIM=c(-0.06,0.06)}else{YLIM=c(0.5,0.7)}

mean_divergence <- t(sapply(data_cell, function(x)   tapply(x$divergence, as.factor(x$class_nb_contact), mean, na.rm=T)))
divergence_conf_low <- t(sapply(data_cell, function(x) tapply(x$divergence, as.factor(x$class_nb_contact), function(y) {z<-t.test(y); return(z[["conf.int"]][1])})))
divergence_conf_high <- t(sapply(data_cell, function(x) tapply(x$divergence, as.factor(x$class_nb_contact), function(y) {z<-t.test(y); return(z[["conf.int"]][2])})))

plot(as.numeric(mean_divergence["Bcell",]), type="l", col=dataset.colors["Bcell"], ylim=YLIM, xlab="", ylab="", axes=F)
lines(as.numeric(mean_divergence["ESC",]), col=dataset.colors["ESC"], lwd=1.5)
lines(as.numeric(mean_divergence["adipo",]), col=dataset.colors["adipo"], lwd=1.5)

## X axis
axis(side=1, at=1:5, labels=1:5, mgp=c(3, 0.65, 0))
axis(side=2, mgp=c(3, 0.75, 0))

## axis labels
mtext("Complexity", side=1, line=2.5, cex=0.9)
mtext("Residual Expression Divergence", side=2, line=2.5, cex=0.9)

## confidence intervals

for(dataset in rownames(mean_divergence)){
  for (col in 1:ncol(mean_divergence)){
    segments(x0=col, x1=col, y0=as.numeric(divergence_conf_low[dataset,col]), y1=as.numeric(divergence_conf_high[dataset,col]), col=dataset.colors[dataset])
  }
}

## legend & plot label
#legend("topright", legend=names(dataset.colors), col=dataset.colors, lty=1, bty='n', inset=c(0.05, 0), xpd=NA, cex=1.4)
mtext("I", side=3, line=1, at=1, font=2, cex=1.2)

########################  G - Expression Divergence vs Number of conserved enhancers ######################## 
mean_divergence <- t(sapply(data_cell, function(x)   tapply(x$divergence, as.factor(x$class_align_score), mean, na.rm=T)))
divergence_conf_low <- t(sapply(data_cell, function(x) tapply(x$divergence, as.factor(x$class_align_score), function(y) {z<-t.test(y); return(z[["conf.int"]][1])})))
divergence_conf_high <- t(sapply(data_cell, function(x) tapply(x$divergence, as.factor(x$class_align_score), function(y) {z<-t.test(y); return(z[["conf.int"]][2])})))

plot(as.numeric(mean_divergence["Bcell",]), type="l", col=dataset.colors["Bcell"], ylim=YLIM, xlab="", ylab="", axes=F)
lines(as.numeric(mean_divergence["ESC",]), col=dataset.colors["ESC"], lwd=1.5)
lines(as.numeric(mean_divergence["adipo",]), col=dataset.colors["adipo"], lwd=1.5)

## X axis
axis(side=1, at=1:5, labels=1:5, mgp=c(3, 0.65, 0))
axis(side=2, mgp=c(3, 0.75, 0))

## axis labels
mtext("Alignment score", side=1, line=2.25, cex=0.9)

## confidence intervals

for(dataset in rownames(mean_divergence)){
  for (col in 1:ncol(mean_divergence)){
    segments(x0=col, x1=col, y0=as.numeric(divergence_conf_low[dataset,col]), y1=as.numeric(divergence_conf_high[dataset,col]), col=dataset.colors[dataset])
  }
}

## plot label
mtext("J", side=3, line=1, at=1, font=2, cex=1)

############## H - Expression Divergence vs Number of conserved contacts ############################# 
mean_divergence <- t(sapply(data_cell, function(x)   tapply(x$divergence, as.factor(x$class_cons_int), mean, na.rm=T)))
divergence_conf_low <- t(sapply(data_cell, function(x) tapply(x$divergence, as.factor(x$class_cons_int), function(y) {z<-t.test(y); return(z[["conf.int"]][1])})))
divergence_conf_high <- t(sapply(data_cell, function(x) tapply(x$divergence, as.factor(x$class_cons_int), function(y) {z<-t.test(y); return(z[["conf.int"]][2])})))

plot(as.numeric(mean_divergence["Bcell",]), type="l", col=dataset.colors["Bcell"], ylim=YLIM, xlab="", ylab="", axes=F)
lines(as.numeric(mean_divergence["ESC",]), col=dataset.colors["ESC"], lwd=1.5)
lines(as.numeric(mean_divergence["adipo",]), col=dataset.colors["adipo"], lwd=1.5)

## X axis
breaks_names=c(">0.1", "0.1-25", "26-50", "51-75", ">75")
axis(side=1, at=1:5, labels=breaks_names, mgp=c(3, 0.65, 0))
axis(side=2, mgp=c(3, 0.75, 0))

## axis labels
mtext("Contacts Conservation", side=1, line=2.25, cex=0.9)

## confidence intervals

for(dataset in rownames(mean_divergence)){
  for (col in 1:ncol(mean_divergence)){
    segments(x0=col, x1=col, y0=as.numeric(divergence_conf_low[dataset,col]), y1=as.numeric(divergence_conf_high[dataset,col]), col=dataset.colors[dataset])
  }
}

## plot label
mtext("K", side=3, line=1, at=1, font=2, cex=1)

#dev.off()
