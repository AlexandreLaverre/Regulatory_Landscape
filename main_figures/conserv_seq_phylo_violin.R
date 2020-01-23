library(ape)
par(mai = c(0.3, 0.2, 0.3, 0.2))
layout(matrix(c(1,1,2,3),nrow = 1,byrow = TRUE))

tree <- read.tree("/home/laverre/Documents/Regulatory_Landscape/data/ensembl_tree")
tree <- keep.tip(tree, c("Mus_musculus", "Homo_sapiens", "Rattus_norvegicus", "Macaca_mulatta", "Oryctolagus_cuniculus", "Canis_lupus_familiaris", "Bos_taurus", "Loxodonta_africana", "Monodelphis_domestica", "Gallus_gallus"))
plot(tree, cex=1.2, y.lim=c(0.4,10.3), x.lim=c(0,1.07), label.offset = 0.01)


load("Fig4_human_mean.Rdata")

vioplot(c(0,conserv_simul), at=c(0,4,8,12,16,20,24,28,32,36), col='cornflowerblue', axes=F, yaxt='n',names=c("human", colnames(conserv)), horizontal = T, las=1, cex.names = 0.8, main="Human conservation")
vioplot(c(0,conserv), at=c(1,5,9,13,17,21,25,29,33,37), col='indianred', add=T, axes=F, horizontal = T)
vioplot(c(0,conserv_enh), at=c(2,6,10,14,18,22,26,30,34,38), col='lightgreen', add=T, axes=F, horizontal = T)
axis(1, at=seq(0,1,0.2), labels=c("0.0", "0.2", "0.4", "0.6", "0.8", "1.0"))

load("Fig4_mouse_mean.Rdata")

vioplot(c(conserv_simul[,c(3,4)],0,conserv_simul[,c(1,2,5,6,7,8,9)]), at=c(0,4,8,12,16,20,24,28,32,36), col='cornflowerblue', axes=F, yaxt='n', horizontal = T, las=1, cex.names = 0.8, main="Mouse conservation")
vioplot(c(conserv[,c(3,4)],0,conserv[,c(1,2,5,6,7,8,9)]), at=c(1,5,9,13,17,21,25,29,33,37), col='indianred', add=T, axes=F, horizontal = T)
vioplot(c(conserv_enh[,c(3,4)],0,conserv_enh[,c(1,2,5,6,7,8,9)]), at=c(2,6,10,14,18,22,26,30,34,38), col='lightgreen', add=T, axes=F, horizontal = T)

axis(1, at=seq(0,1,0.2), labels=c("0.0", "0.2", "0.4", "0.6", "0.8", "1.0"))
