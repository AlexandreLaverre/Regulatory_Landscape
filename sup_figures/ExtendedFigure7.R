######################################################################################################################

library(Hmisc)

source("../main_figures/parameters.R")

######################################################################################################################

sp="human"

load(paste(pathFigures, "RData/", sp, ".cells.types.parallel.trends.RData", sep=""))

cells <- c("Bcell", "ESC", "adipo")
dataset.colors=c("firebrick1", "forestgreen", "navy")
names(dataset.colors) = cells

pdf(file=paste(pathFigures, "/ExtendedFigure7.pdf", sep=""), width=7, height=6)

m=matrix(rep(NA, 2*13), nrow=2)
m[1,]=c(rep(1,4),  rep(2,3),  rep(3,3), rep(4,3))
m[2,]=c(rep(5,4), rep(6,3),  rep(7,3), rep(8,3))
layout(m)

######################################################################################################################
####################################### Parallel trends among cell types #############################################
par(mai = c(0.7, 0.1, 0.5, 0.1)) # bottom, left, top, right
cex.label=1

##  Correlation of expression levels

dotchart(correl_expression[rev(cells), "Spearman"], col=dataset.colors[rev(cells)], labels=c("pre-\nadipocytes", "ESC", "B cells"), pch=20, pt.cex=1.1,  xlim=c(min(correl_expression)-0.02, max(correl_expression)+0.02), cex.lab=1.25)
mtext("expression correlation\n(Spearman's rho)", side=1, line=3.5, cex=0.7)

## plot label
mtext("a", side=3, at=0.71, font=2, cex=cex.label, line=1.5)

#  Conserved Expression
dotchart(conserv_expression[rev(cells),"Mean"], col=dataset.colors[rev(cells)], pch=20, labels="", pt.cex=1.1,
         xlim=c(min(conserv_expression[,"Conf_low"])-0.01, max(conserv_expression[,"Conf_high"])+0.01))
segments(x0=conserv_expression[rev(cells),"Conf_low"], x1=conserv_expression[rev(cells),"Conf_high"], y0=1:3, y1=1:3, col=dataset.colors[rev(cells)], lwd=1.1)
## segments(x0=conserv_expression[rev(cells),"Conf_low"], x1=conserv_expression[rev(cells),"Conf_low"], y0=(1:3)-0.02, y1=(1:3)+0.02, col=dataset.colors[rev(cells)], lwd=1)
## segments(x0=conserv_expression[rev(cells),"Conf_high"], x1=conserv_expression[rev(cells),"Conf_high"], y0=(1:3)-0.02, y1=(1:3)+0.02, col=dataset.colors[rev(cells)], lwd=1)

mtext("expression conservation\n(corrected)", side=1, line=3.5, cex=0.7)

## plot label
mtext("b", side=3, at=-0.11, font=2, cex=cex.label, line=1.5)

#  Correlattion of complexity zscore 
dotchart(correl_complexity[rev(cells),"Spearman"], col=dataset.colors[rev(cells)], labels='', pch=20, pt.cex=1.1, 
         xlim=c(min(correl_complexity[,"Spearman"])-0.01, max(correl_complexity[,"Spearman"])+0.01))
mtext("conservation of regulatory complexity\n(Spearman's rho)", side=1, line=3.5, cex=0.7)

## plot label
mtext("c", side=3, at=0.215, font=2, cex=cex.label, line=1.5)

# #  dN / dS
dotchart(gene_dnds[rev(cells),"Mean"], col=dataset.colors[rev(cells)], labels='', pch=20, pt.cex=1.1,
         xlim=c(min(gene_dnds[,"Conf_low"])-0.02, max(gene_dnds[,"Conf_high"])+0.02))
segments(x0=gene_dnds[rev(cells),"Conf_low"], x1=gene_dnds[rev(cells),"Conf_high"], y0=1:3, y1=1:3, col=dataset.colors[rev(cells)], lwd=1.1)
## segments(x0=gene_dnds[rev(cells),"Conf_low"], x1=gene_dnds[rev(cells),"Conf_low"], y0=(1:3)-0.02, y1=(1:3)+0.02, col=dataset.colors[rev(cells)], lwd=1)
## segments(x0=gene_dnds[rev(cells),"Conf_high"], x1=gene_dnds[rev(cells),"Conf_high"], y0=(1:3)-0.02, y1=(1:3)+0.02, col=dataset.colors[rev(cells)], lwd=1)

mtext("1- dN/dS \nof top expressed genes", side=1, line=3.5, cex=0.7)

## plot label
mtext("d", side=3, at=0.855, font=2, cex=cex.label, line=1.5)

#  Enhancer Alignment score
dotchart(enh_evol[rev(cells),"Mean"], col=dataset.colors[rev(cells)], pch=20, labels=c("pre-\nadipocytes", "ESC", "B cells"), pt.cex=1.1,
         xlim=c(min(enh_evol[,"Conf_low"])-0.02, max(enh_evol[,"Conf_high"])+0.02))
segments(x0=enh_evol[rev(cells),"Conf_low"], x1=enh_evol[rev(cells),"Conf_high"], y0=1:3, y1=1:3, col=dataset.colors[rev(cells)])
## segments(x0=enh_evol[rev(cells),"Conf_low"], x1=enh_evol[rev(cells),"Conf_low"], y0=(1:3)-0.02, y1=(1:3)+0.02, col=dataset.colors[rev(cells)], lwd=1)
## segments(x0=enh_evol[rev(cells),"Conf_high"], x1=enh_evol[rev(cells),"Conf_high"], y0=(1:3)-0.02, y1=(1:3)+0.02, col=dataset.colors[rev(cells)], lwd=1)

mtext("alignment score\n of contacted enhancers", side=1, line=3.5, cex=0.7)

## plot label
mtext("e", side=3, at=0.4725, font=2, cex=cex.label, line=1.5)

# Conserved sequence
dotchart(seq_conserv[rev(cells),"Mean"], col=dataset.colors[rev(cells)], pch=20, labels='', pt.cex=1.1,
         xlim=c(min(seq_conserv[,"Conf_low"])-0.02, max(seq_conserv[,"Conf_high"])+0.02))
segments(x0=seq_conserv[rev(cells),"Conf_low"], x1=seq_conserv[rev(cells),"Conf_high"], y0=1:3, y1=1:3, col=dataset.colors[rev(cells)])
## segments(x0=seq_conserv[rev(cells),"Conf_low"], x1=seq_conserv[rev(cells),"Conf_low"], y0=(1:3)-0.02, y1=(1:3)+0.02, col=dataset.colors[rev(cells)], lwd=1)
## segments(x0=seq_conserv[rev(cells),"Conf_high"], x1=seq_conserv[rev(cells),"Conf_high"], y0=(1:3)-0.02, y1=(1:3)+0.02, col=dataset.colors[rev(cells)], lwd=1)

mtext("proportion of conserved\n enhancers by gene", side=1, line=3.5, cex=0.7)

## plot label
mtext("f", side=3, at=0.53, font=2, cex=cex.label, line=1.5)

#  Conserved synteny
dotchart(synteny_conserv[rev(cells),"Mean"], col=dataset.colors[rev(cells)], pch=20, labels='', pt.cex=1.1,
         xlim=c(min(synteny_conserv[,"Conf_low"])-0.04, max(synteny_conserv[,"Conf_high"])+0.04))
segments(x0=synteny_conserv[rev(cells),"Conf_low"], x1=synteny_conserv[rev(cells),"Conf_high"], y0=1:3, y1=1:3, col=dataset.colors[rev(cells)], lwd=1)
## segments(x0=synteny_conserv[rev(cells),"Conf_low"], x1=synteny_conserv[rev(cells),"Conf_low"], y0=(1:3)-0.02, y1=(1:3)+0.02, col=dataset.colors[rev(cells)], lwd=1)
## segments(x0=synteny_conserv[rev(cells),"Conf_high"], x1=synteny_conserv[rev(cells),"Conf_high"], y0=(1:3)-0.02, y1=(1:3)+0.02, col=dataset.colors[rev(cells)], lwd=1)

mtext("proportion of enhancers\nmaintained in synteny", side=1, line=3.5, cex=0.7)
## plot label
mtext("g", side=3, at=0.85, font=2, cex=cex.label, line=1.5)

#  Conserved contacts
dotchart(contact_conserv[rev(cells),"Mean"], col=dataset.colors[rev(cells)], pch=20, labels='', pt.cex=1.1,
         xlim=c(min(contact_conserv[,"Conf_low"])-0.04, max(contact_conserv[,"Conf_high"])+0.04))
segments(x0=contact_conserv[rev(cells),"Conf_low"], x1=contact_conserv[rev(cells),"Conf_high"], y0=1:3, y1=1:3, col=dataset.colors[rev(cells)], lwd=1)
## segments(x0=contact_conserv[rev(cells),"Conf_low"], x1=contact_conserv[rev(cells),"Conf_low"], y0=(1:3)-0.02, y1=(1:3)+0.02, col=dataset.colors[rev(cells)], lwd=1)
## segments(x0=contact_conserv[rev(cells),"Conf_high"], x1=contact_conserv[rev(cells),"Conf_high"], y0=(1:3)-0.02, y1=(1:3)+0.02, col=dataset.colors[rev(cells)], lwd=1)

mtext("percentage of enhancers\nmaintained in contact", side=1, line=3.5, cex=0.7)

## plot label
mtext("h", side=3, at=0.12, font=2, cex=cex.label, line=1.5)

##############################################################################################################

dev.off()

##############################################################################################################
