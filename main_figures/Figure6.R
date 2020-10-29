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
#pdf(file=paste(pathFigures, "Figure6.pdf", sep=""))

# m=matrix(rep(NA, 2*11), nrow=2)
# m[1,]=c(rep(1,2), rep(2, 3), rep(3,3), rep(4,3))
# m[2,]=c(rep(5,2), rep(6, 3), rep(7,3), rep(8,3))
# layout(m)
# par(mai = c(0, 0, 0, 0)) # bottom, left, top, right
# plot.new()
# legend("center", legend=label.enhancers[enhancer.datasets[[sp]]], pch=20, 
#        col=col.enhancers[enhancer.datasets[[sp]]], cex=1.1, y.intersp = 1.5,
#        bty="o", box.col="white", bg="white",  inset=c(0, 0.03))


par(mfrow=c(2,3))

par(mai = c(0.5, 0.5, 0.5, 0)) # bottom, left, top, right

################################################################################################################################
############################## PART1 : All cells & Cardoso-Moreira  ##########################################################
if (DivergenceMeasure == "CorrelationSpearman"){ylim=c(0.45, 0.65)}else{ylim=c(-0.03, 0.1)}

#### A - Gene expression profil similarity and enhancers conserved in sequences ####
xpos=seq(1, length(levels(expdiv_all$classTau)), 1)
names(xpos) = levels(expdiv_all$classTau)

smallx=c(-0.15, -0.075, 0.075, 0.15)
names(smallx)=enhancer.datasets[[sp]]

xlim=c(0.5, length(levels(regland[[enh]]$class_align_score))+0.5)

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

axis(side=2, mgp=c(3, 0.75, 0), cex.axis=1.1)
mtext("Spearman's coefficient", side=2, line=2.5, cex=0.9)

mtext("A", side=3, at=0.45, font=2, cex=1.1, line=0.5)
legend("bottomright", legend=label.enhancers[enhancer.datasets[[sp]]], pch=20, 
             col=col.enhancers[enhancer.datasets[[sp]]], cex=1,
             bty="o", box.col="white", bg="white",  inset=c(0, 0.01))

#### B - Gene expression profil similarity and enhancers conserved in synteny ####
par(mai = c(0.5, 0.2, 0.5, 0)) # bottom, left, top, right

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
mtext("B", side=3, at=0.45, font=2, cex=1.1, line=0.5)

#### C - Gene expression profil similarity and enhancers conserved in contact ####
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
mtext("C", side=3, at=0.45, font=2, cex=1.1, line=0.5)


################################################################################################################################
######################################### PART2 : Common cells types  ##########################################################
# par(mai = c(0, 0, 0, 0)) # bottom, left, top, right
# plot.new()
# legend("center", legend=c("Bcell", "ESC", "Pre-adipocytes"),
#        lty=1, col=dataset.colors[cells], bty="n", cex=1.2, y.intersp = 2)

par(mai = c(0.5, 0.5, 0.5, 0)) # bottom, left, top, right
YLIM=c(0.3,0.6)

########################  D - Expression Conservation vs Number of conserved enhancers ######################## 
mean_divergence <- t(sapply(data_cell, function(x)   tapply(x$Conservation, as.factor(x$class_align_score), mean, na.rm=T)))
divergence_conf_low <- t(sapply(data_cell, function(x) tapply(x$Conservation, as.factor(x$class_align_score), function(y) {z<-t.test(y); return(z[["conf.int"]][1])})))
divergence_conf_high <- t(sapply(data_cell, function(x) tapply(x$Conservation, as.factor(x$class_align_score), function(y) {z<-t.test(y); return(z[["conf.int"]][2])})))

plot(as.numeric(mean_divergence["Bcell",]), type="l", col=dataset.colors["Bcell"], ylim=YLIM, xlab="", ylab="", axes=F)
lines(as.numeric(mean_divergence["ESC",]), col=dataset.colors["ESC"], lwd=1.5)
lines(as.numeric(mean_divergence["adipo",]), col=dataset.colors["adipo"], lwd=1.5)

## X axis
axis(side=1, at=1:5, labels=1:5, mgp=c(3, 0.65, 0))
axis(side=2, mgp=c(3, 0.75, 0))

## axis labels
mtext("Alignment score", side=1, line=2.25, cex=0.9)
mtext("Expression Conservation", side=2, line=2.5, cex=0.9)

#legend
legend("topright", legend=c("Bcell", "ESC", "Pre-adipocytes"), lty=1, col=dataset.colors[cells], bty="n", cex=1.2)

## confidence intervals
for(dataset in rownames(mean_divergence)){
  for (col in 1:ncol(mean_divergence)){
    segments(x0=col, x1=col, y0=as.numeric(divergence_conf_low[dataset,col]), y1=as.numeric(divergence_conf_high[dataset,col]), col=dataset.colors[dataset])
  }
}

## plot label
mtext("D", side=3, line=1, at=1, font=2, cex=1)

############## E - Expression Conservation vs Number of conserved contacts ############################# 
mean_divergence <- t(sapply(data_cell, function(x)   tapply(x$Conservation, as.factor(x$class_cons_synt), mean, na.rm=T)))
divergence_conf_low <- t(sapply(data_cell, function(x) tapply(x$Conservation, as.factor(x$class_cons_synt), function(y) {z<-t.test(y); return(z[["conf.int"]][1])})))
divergence_conf_high <- t(sapply(data_cell, function(x) tapply(x$Conservation, as.factor(x$class_cons_synt), function(y) {z<-t.test(y); return(z[["conf.int"]][2])})))

plot(as.numeric(mean_divergence["Bcell",]), type="l", col=dataset.colors["Bcell"], ylim=YLIM, xlab="", ylab="", axes=F)
lines(as.numeric(mean_divergence["ESC",]), col=dataset.colors["ESC"], lwd=1.5)
lines(as.numeric(mean_divergence["adipo",]), col=dataset.colors["adipo"], lwd=1.5)

## X axis
breaks_names=c(">75%", "75-99%", ">99%")
axis(side=1, at=1:3, labels=breaks_names, mgp=c(3, 0.65, 0))

## axis labels
mtext("Synteny Conservation", side=1, line=2.25, cex=0.9)

## confidence intervals

for(dataset in rownames(mean_divergence)){
  for (col in 1:ncol(mean_divergence)){
    segments(x0=col, x1=col, y0=as.numeric(divergence_conf_low[dataset,col]), y1=as.numeric(divergence_conf_high[dataset,col]), col=dataset.colors[dataset])
  }
}

## plot label
mtext("E", side=3, line=1, at=1, font=2, cex=1)

############## F - Expression Conservation vs Number of conserved contacts ############################# 
mean_divergence <- t(sapply(data_cell, function(x)   tapply(x$Conservation, as.factor(x$class_cons_int), mean, na.rm=T)))
divergence_conf_low <- t(sapply(data_cell, function(x) tapply(x$Conservation, as.factor(x$class_cons_int), function(y) {z<-t.test(y); return(z[["conf.int"]][1])})))
divergence_conf_high <- t(sapply(data_cell, function(x) tapply(x$Conservation, as.factor(x$class_cons_int), function(y) {z<-t.test(y); return(z[["conf.int"]][2])})))

plot(as.numeric(mean_divergence["Bcell",]), type="l", col=dataset.colors["Bcell"], ylim=YLIM, xlab="", ylab="", axes=F)
lines(as.numeric(mean_divergence["ESC",]), col=dataset.colors["ESC"], lwd=1.5)
lines(as.numeric(mean_divergence["adipo",]), col=dataset.colors["adipo"], lwd=1.5)

## X axis
breaks_names=c("<1%", "1-25%", "25-50%", "50-75%", ">75%")
axis(side=1, at=1:5, labels=breaks_names, mgp=c(3, 0.65, 0))

## axis labels
mtext("Contacts Conservation", side=1, line=2.25, cex=0.9)

## confidence intervals

for(dataset in rownames(mean_divergence)){
  for (col in 1:ncol(mean_divergence)){
    segments(x0=col, x1=col, y0=as.numeric(divergence_conf_low[dataset,col]), y1=as.numeric(divergence_conf_high[dataset,col]), col=dataset.colors[dataset])
  }
}

## plot label
mtext("F", side=3, line=1, at=1, font=2, cex=1)

#dev.off()
