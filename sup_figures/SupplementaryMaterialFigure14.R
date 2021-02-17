###########################################################################

objects=ls()

if(!"pathFigures"%in%objects){

 source("../main_figures/parameters.R")

 library(ape)
 library(vioplot)

 load=T
 prepare=T
}

###########################################################################


if(load){
 ref_sp = "mouse"
 
 enhancers = enhancer.datasets[[ref_sp]]

 load(paste(pathFigures, "RData/data.sequence.conservation.stats.pcungapped.", ref_sp, ".RData", sep=""))

 load=F
}

#########################################################################################################################

## 1 column width 85 mm = 3.34 in
## 1.5 column width 114 mm = 4.49 in
## 2 columns width 174 mm = 6.85 in
## max height: 11 in

#########################################################################################################################

pdf(paste(pathFigures, "SupplementaryMaterialFigure13.pdf", sep=""), width=6.85, height=4)

par(mai = c(0.5, 0.1, 0.3, 0.1)) #bottom, left, top and right

m=matrix(rep(NA, 3*8), nrow=3)

for(i in 1:3){
  m[i,]=c(1,1,2,2,3,3, 4, 4)
}

layout(m)

################################## a - Phylogenetic tree ################################################################

tree <- read.tree(paste(pathFigures, "RData/Ensembl_species_tree", sep=""))

tree <- keep.tip(tree, c("Mus_musculus", "Homo_sapiens", "Rattus_norvegicus", "Macaca_mulatta", "Oryctolagus_cuniculus", "Canis_lupus_familiaris", "Bos_taurus", "Loxodonta_africana", "Monodelphis_domestica", "Gallus_gallus"))

species_names <- c("human", "macaque", "mouse", "rat", "rabbit", "cow", "dog", "elephant", "opossum", "chicken")
names(species_names) <- c("Homo_sapiens", "Macaca_mulatta", "Mus_musculus", "Rattus_norvegicus",  "Oryctolagus_cuniculus",  "Bos_taurus", "Canis_lupus_familiaris","Loxodonta_africana", "Monodelphis_domestica", "Gallus_gallus")
  
par(mar=c(4.1,1.1, 2.1, 1))
plot(tree, cex=1.1, y.lim=c(0.5,10.5), x.lim=c(0,1.07), label.offset = 0.01, show.tip.label = F, main="")
tiplabels(species_names[tree$tip.label], bg = NA, adj = -0.1, frame="none", cex=1.1, xpd=NA)

 # legend for the plot
legend("bottomleft", fill=dataset.colors, border=dataset.colors, legend = c("PCHi-C data", "simulated data"), bty='n', cex=1.1, xpd=T, inset=c(-0.01, -0.15), horiz=FALSE)

# label
mtext("a", side=3, line=0.5, at=-0.05, font=2, cex=1.05)

######################## b - Restriction fragments sequence conservation ########################
species <-c("human", "macaque", "rat", "rabbit", "cow", "dog", "elephant", "opossum", "chicken")

ylim=c(0, 40.5)
xlim=c(0, 100)

ypos.sim=2+c(0, 4, 12,16,20,24,28,32,36)-0.27
ypos.obs=2+c(1, 5, 13,17,21,25,29,33,37)+0.27

plot(1, type="n", xlab="", ylab="", axes=F, ylim=ylim, xlim=xlim, main="", bty="n")

# Simulated
par(bty="n")
vioplot(100*frag_align_simul[,species], at=ypos.sim, add=T, border=dataset.colors["Simulated"], col=rgb(t(col2rgb(dataset.colors["Simulated"])/255), alpha = 0.6), plotCentre="line", axes=F, xaxt="n", yaxt="n", horizontal = T, las=1, cex.main = 1.2, main="", bty="n")

# Original
vioplot(100*frag_align_obs[,species], at=ypos.obs, add=T, axes=F, xaxt="n", yaxt="n", horizontal = T, border=dataset.colors["Original"], col=rgb(t(col2rgb(dataset.colors["Original"])/255), alpha = 0.6), plotCentre="line")

## add mean point
points(x=100*apply(frag_align_simul[,species], 2, mean, na.rm=T), y = ypos.sim, col = "white", pch=20, cex=0.8)
points(x=100*apply(frag_align_obs[,species], 2, mean, na.rm=T), y = ypos.obs, col = "white", pch=20, cex=0.8)

## axis and legend
axis(1, pos=0.7, at=seq(0,100,20), labels=c("0", "20", "40", "60", "80", "100"), cex.axis=1, mgp=c(3, 0.75, 0))
mtext("% aligned sequence", side=1, xpd = TRUE, cex=0.75, line=1.25)
mtext("restriction fragments", side=3, line=-0.5, cex=0.75)
mtext("b", side=3, line=0.5, at=-8, font=2, cex=1.05)

######################## c - enhancer sequence conservation ########################

labels=letters[3:4]
names(labels)=c("ENCODE", "FANTOM5")

for(enh in c("ENCODE", "FANTOM5")){
  align_enhancer_obs=list_align_enh[[enh]][["enh_align_obs"]]
  align_enhancer_simul=list_align_enh[[enh]][["enh_align_simul"]]

  plot(1, type="n", xlab="", ylab="", axes=F, ylim=ylim, xlim=xlim, main="", bty="n")
  
   # simulated enhancers
  vioplot(100*align_enhancer_simul[,species], at=ypos.sim, col=rgb(t(col2rgb(dataset.colors["Simulated"])/255), alpha = 0.6), border=dataset.colors["Simulated"], add=T, axes=F, xaxt="n", yaxt="n", horizontal = T, las=1, cex.main = 1.2, main="", plotCentre="line")

   # original enhancer
  vioplot(100*align_enhancer_obs[,species],at=ypos.obs, col=rgb(t(col2rgb(dataset.colors["Original"])/255), alpha = 0.6), border=dataset.colors["Original"], add=T, axes=F, xaxt="n", yaxt="n", horizontal = T, plotCentre="line")
  
# add mean point
  points(x = apply(100*align_enhancer_simul[,species], 2, mean, na.rm=T), y=ypos.sim, col = "white", pch=20, cex=0.8)
  points(x = apply(100*align_enhancer_obs[,species], 2, mean, na.rm=T), y = ypos.obs, col = "white", pch=20, cex=0.8)
  
  # axis and legend
  axis(1, pos=0.7, at=seq(0,100,20), labels=c("0", "20", "40", "60", "80", "100"), cex.axis=1, mgp=c(3, 0.75, 0))
  
  mtext("% aligned sequence", side=1, xpd = TRUE, cex=0.75, line=1.25)
  mtext(paste(enh.syn[enh], "enhancers"), side=3, line=-0.5, cex=0.75)

  ## plot label
  mtext(labels[enh], side=3, line=0.5, at=-8, font=2, cex=1.05)

}

dev.off()

#######################################################################################################
