######################################################################################################################
library(Hmisc)

#setwd("/home/laverre/Data/Regulatory_landscape/scripts/main_figures")
source("parameters.R") ## pathFinalData are defined based on the user name

sp="human"

load(paste(pathFigures, "RData/Fig6_", sp, "_common_cells.Rdata", sep=""))

cells <- c("Bcell", "ESC", "adipo")
dataset.colors=c("firebrick1", "forestgreen", "navy")
names(dataset.colors) = cells

pdf(file=paste(pathFigures, "/SupplementaryFigureX_ParallelTrendsCellTypes.pdf", sep=""), width=7, height=3)
#par(mfrow=c(1,6))
m=matrix(rep(NA, 11), nrow=1)
m[1,]=c(1, rep(2,2),  rep(3,2),  rep(4,2), rep(5,2), rep(6,2))
layout(m)
par(mai = c(0.8, 0, 0.8, 0)) # bottom, left, top, right

######################################################################################################################
####################################### Parallel trends among cell types #############################################
#  Legend
plot.new()
legend("center",col=dataset.colors[cells], legend = cells, bty='n', pch=16, y.intersp=4, cex=1.2)

par(mai = c(0.8, 0.1, 0.5, 0.1)) # bottom, left, top, right
#  Correlation of expression level
dotchart(correl_expression[rev(cells),"Pearson"], col=dataset.colors[rev(cells)], labels='', pch=16, pt.cex=0.8,
         xlim=c(min(correl_expression)-0.02, max(correl_expression)+0.02))
mtext("Expression level \n correlation (rho)", side=1, line=3.5, cex=0.8)
mtext("A", side=3, at=0.7, font=2, cex=1.1, line=1.5)

#  Correlattion of complexity zscore 
dotchart(correl_complexity[rev(cells),"Pearson"], col=dataset.colors[rev(cells)], labels='', pch=16, pt.cex=0.8, 
         xlim=c(min(correl_complexity), max(correl_complexity)+0.01))
mtext("Complexity \n correlation (rho)", side=1, line=3.5, cex=0.8)
mtext("B", side=3, at=0.22, font=2, cex=1.1, line=1.5)

#  dN / dS
dotchart(gene_dnds[rev(cells),"Mean"], col=dataset.colors[rev(cells)], labels='', pch=16, pt.cex=0.8, 
         xlim=c(min(gene_dnds[,"Conf_low"])-0.02, max(gene_dnds[,"Conf_high"])+0.02))
segments(x0=gene_dnds[rev(cells),"Conf_low"], x1=gene_dnds[rev(cells),"Conf_high"], y0=1:3, y1=1:3, col=dataset.colors[rev(cells)])
segments(x0=gene_dnds[rev(cells),"Conf_low"], x1=gene_dnds[rev(cells),"Conf_low"], y0=(1:3)-0.02, y1=(1:3)+0.02, col=dataset.colors[rev(cells)], lwd=1)
segments(x0=gene_dnds[rev(cells),"Conf_high"], x1=gene_dnds[rev(cells),"Conf_high"], y0=(1:3)-0.02, y1=(1:3)+0.02, col=dataset.colors[rev(cells)], lwd=1)

mtext("dN/dS", side=1, line=2.5, cex=0.8)
mtext("C", side=3, at=0.075, font=2, cex=1.1, line=1.5)

#  Enhancer Alignment
dotchart(enh_evol[rev(cells),"Mean"], col=dataset.colors[rev(cells)], pch=16, labels='', pt.cex=0.8,
         xlim=c(min(enh_evol[,"Conf_low"])-0.02, max(enh_evol[,"Conf_high"])+0.02))
segments(x0=enh_evol[rev(cells),"Conf_low"], x1=enh_evol[rev(cells),"Conf_high"], y0=1:3, y1=1:3, col=dataset.colors[rev(cells)])
segments(x0=enh_evol[rev(cells),"Conf_low"], x1=enh_evol[rev(cells),"Conf_low"], y0=(1:3)-0.02, y1=(1:3)+0.02, col=dataset.colors[rev(cells)], lwd=1)
segments(x0=enh_evol[rev(cells),"Conf_high"], x1=enh_evol[rev(cells),"Conf_high"], y0=(1:3)-0.02, y1=(1:3)+0.02, col=dataset.colors[rev(cells)], lwd=1)

mtext("Alignment Score", side=1, line=2.5, cex=0.8)
mtext("D", side=3, at=0.4, font=2, cex=1.1, line=1.5)

#  Conserved contacts
dotchart(contact_conserv[rev(cells),"Mean"], col=dataset.colors[rev(cells)], pch=16, labels='', pt.cex=0.8,
         xlim=c(min(contact_conserv[,"Conf_low"])-0.04, max(contact_conserv[,"Conf_high"])+0.04))
segments(x0=contact_conserv[rev(cells),"Conf_low"], x1=contact_conserv[rev(cells),"Conf_high"], y0=1:3, y1=1:3, col=dataset.colors[rev(cells)], lwd=1)
segments(x0=contact_conserv[rev(cells),"Conf_low"], x1=contact_conserv[rev(cells),"Conf_low"], y0=(1:3)-0.02, y1=(1:3)+0.02, col=dataset.colors[rev(cells)], lwd=1)
segments(x0=contact_conserv[rev(cells),"Conf_high"], x1=contact_conserv[rev(cells),"Conf_high"], y0=(1:3)-0.02, y1=(1:3)+0.02, col=dataset.colors[rev(cells)], lwd=1)

mtext("Conserved\ncontact", side=1, line=3.5, cex=0.8)
mtext("E", side=3, at=0.03, font=2, cex=1.1, line=1.5)
#########
dev.off()
