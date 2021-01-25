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
 ref_sp = "human"
 
 enhancers = enhancer.datasets[[ref_sp]]

 load(paste(pathFigures, "RData/data.sequence.conservation.pcungapped.", ref_sp, ".Rdata", sep=""))

 load=F
}

#########################################################################################################################

## 1 column width 85 mm = 3.34 in
## 1.5 column width 114 mm = 4.49 in
## 2 columns width 174 mm = 6.85 in
## max height: 11 in

#########################################################################################################################

pdf(paste(pathFigures, "SupplementaryFigure13.pdf", sep=""), width=6.85, height=4)

par(mai = c(0.5, 0.1, 0.3, 0.1)) #bottom, left, top and right

m=matrix(rep(NA, 3*8), nrow=3)

for(i in 1:3){
  m[i,]=c(1,1,2,2,3,3, 4, 4)
}

layout(m)

################################## a - Phylogenetic tree ################################################################

tree <- read.tree(paste(pathFigures, "RData/Ensembl_species_tree", sep=""))

tree <- keep.tip(tree, c("Mus_musculus", "Homo_sapiens", "Rattus_norvegicus", "Macaca_mulatta", "Oryctolagus_cuniculus", "Canis_lupus_familiaris", "Bos_taurus", "Loxodonta_africana", "Monodelphis_domestica", "Gallus_gallus"))
species <- c("macaque", "mouse", "rat", "rabbit", "cow", "dog", "elephant", "opossum", "chicken")
species_names <- c("human", species)

par(mar=c(4.1,1.1, 2.1, 1))
plot(tree, cex=1.1, y.lim=c(0.3,10.5), x.lim=c(0,1.07), label.offset = 0.01, show.tip.label = F, main="")
tiplabels(species_names, bg = NA, adj = -0.1, frame="none", cex=1.1, xpd=NA)

 # legend for the plot

legend("bottomleft", fill=dataset.colors, border=dataset.colors, legend = c("PCHi-C data", "simulated data"), bty='n', cex=1.2, xpd=T, inset=c(-0.01, -0.15), horiz=FALSE)

# label
mtext("a", side=3, line=0.5, at=-0.05, font=2, cex=1.2)

######################## b - enhancer sequence conservation ########################


ylim=c(-2, 38.5)
xlim=c(0, 100)
ypos.sim=c(4,8,12,16,20,24,28,32,36)-0.27
ypos.obs=c(5,9,13,17,21,25,29,33,37)+0.27


labels=letters[2:4]
names(labels)=c("FANTOM5", "RoadmapEpigenomics", "FOCS_GRO_seq")

for(enh in c("FANTOM5", "RoadmapEpigenomics", "FOCS_GRO_seq")){
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
  axis(1, pos=0.7, at=seq(0,100,20), labels=c("0", "20", "40", "60", "80", "100"), cex.axis=1)
  
  mtext("% aligned sequence", side=1, xpd = TRUE, cex=0.75, line=0.5)
  
  mtext(enh.syn[enh], side=3, line=-1, cex=0.7)
  
  mtext(labels[enh], side=3, line=0.5, at=-8, font=2, cex=1.2)

}

dev.off()

#######################################################################################################
