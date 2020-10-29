######################################################################################################################
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

######################################################################################################################
##### A - Regulatory Landscape Complexity and average gene expression level across samples ####

##### B - Regulatory Landscape Complexity and  gene expression level across common cell types ####

##### C - Regulatory Landscape Complexity and  gene expression specificity (Cardoso-Moreira) ####
xpos=seq(1, length(levels(expdiv_all$classTau)), 1)
names(xpos) = levels(expdiv_all$classTau)

smallx=c(-0.15, -0.075, 0.075, 0.15)
names(smallx)=enhancer.datasets[[sp]]

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

##### D - Regulatory Landscape Complexity and  gene expression profil similarity (Spearman) ####
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

##### E - Regulatory Landscape Complexity and  gene expression level evolution in common cell types ####
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
