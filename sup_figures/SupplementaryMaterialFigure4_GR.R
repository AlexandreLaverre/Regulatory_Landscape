#################################################################################################################
## if it's the first time we run this figure, we load and prepare data

objects=ls()

if(!"pathScripts"%in%objects){
  load=T
  prepare=T
  source("../main_figures/parameters.R")
  
 
}


#################################################################################################################

if(load){
  ref_sp="mouse"
   
  load(paste(pathFigures, "RData/data.promoter.enhancer.correlation.",ref_sp,".RData", sep=""))
  
  enhancers = enhancer.datasets[[ref_sp]]

  all.enh.prop=data.frame(id=character(0), data=numeric(0), conf_up=numeric(0), conf_low=numeric(0))
  all.enh.prop.nbcell=list("obs"=list(), "simul"=list())
  all.enh.prop.dist=list("obs"=list(), "simul"=list())
  
  for(enh in enhancers){
    load(paste(pathFigures, "RData/data.enhancer.coverage.", ref_sp,".",enh,".RData", sep=""))

    all.enh.prop=rbind(all.enh.prop, enh_prop, stringsAsFactors=F)
    
    all.enh.prop.nbcell[["obs"]]=c(all.enh.prop.nbcell[["obs"]], enh_prop_nb_cell[["obs"]])
    all.enh.prop.nbcell[["simul"]]=c(all.enh.prop.nbcell[["simul"]], enh_prop_nb_cell[["simul"]])

    all.enh.prop.dist[["obs"]]=c(all.enh.prop.dist[["obs"]], enh_prop_dist[["obs"]])
    all.enh.prop.dist[["simul"]]=c(all.enh.prop.dist[["simul"]], enh_prop_dist[["simul"]])
  }

  enh_prop=all.enh.prop
  enh_prop_dist=all.enh.prop.dist
  enh_prop_nb_cell=all.enh.prop.nbcell
}

#################################################################################################################

## 1 column width 85 mm = 3.34 in
## 1.5 column width 114 mm = 4.49 in 
## 2 columns width 174 mm = 6.85 in
## max height: 11 in

#################################################################################################################

pdf(paste(pathFigures, "GenomeResearch_Figures/SupplementaryMaterialFigure4.pdf", sep=""), width=6.85, height=5)

par(mai = c(0.5, 0.5, 0.3, 0.2)) # bottom, left, top, right

m=matrix(rep(NA, 2*10), nrow=2)
m[1,]=c(rep(1,5), rep(2,5))
m[2,]=c(rep(3,5), rep(4,5))
layout(m)

##################################################################################################################
############################################  % length covered by enhancers  #####################################

m.prop=t(matrix(enh_prop$data, nrow=length(enhancers), byrow=T))

par(mar=c(2.1, 4.5, 2.75, 1))

barcenter <- barplot(m.prop, beside=T,  border=dataset.colors, col=dataset.colors, 
                     lwd=1.5, cex.names=0.8,
                     ylim=c(0,15), ylab="", xlab="", axisnames = F, main="", space=c(0.1, 1), las=2)


xposlab=apply(barcenter, 2, mean)
allxpos=as.numeric(barcenter)

lab = enhancers

mtext(lab, line=c(rep(0.5,length(enhancers)/2), rep(1.3,length(enhancers)/2)), side=1, at=xposlab, cex=0.75)
mtext("% length covered by enhancers", side=2, cex=0.85, line=2.7, at=7)

at="topright"
legend(at, legend = c("PCHi-C data", "simulated data"), fill=dataset.colors, border=dataset.colors,  bty='n', cex=1.1, inset=c(0, -0.1), xpd=NA)

par(lwd=1)
segments(allxpos, enh_prop$conf_up, allxpos, enh_prop$conf_low, lwd = 3)
arrows(allxpos, enh_prop$conf_up, allxpos, enh_prop$conf_low, lwd = 1.5, angle = 90, code = 3, length = 0.05)

for (x in seq(1,length(allxpos)-1, by=2)){
  segments(allxpos[x], enh_prop$data[x]+1, allxpos[x+1], enh_prop$data[x]+1) 
  text("***", x=(allxpos[x]+allxpos[x+1])/2, y=enh_prop$data[x]+1.5, cex=1.2)
}

at=-0.1
mtext("A", side=3, line=1.45, at=at, font=2, cex=1.2)

###########################################################################################################################
##########################  Fig2-B - Enhancer proportion according to distance ############################################

## only one enhancer dataset

enh="FANTOM5" ## for this dataset we have promoter-enhancer correlations

ymax=max(c(enh_prop_dist[["obs"]][[paste0(enh,"_conflow")]], enh_prop_dist[["obs"]][[paste0(enh,"_confup")]],  enh_prop_dist[["simul"]][[paste0(enh,"_conflow")]], enh_prop_dist[["simul"]][[paste0(enh,"_confup")]]))
ymin=min(c(enh_prop_dist[["obs"]][[paste0(enh,"_conflow")]], enh_prop_dist[["obs"]][[paste0(enh,"_confup")]],  enh_prop_dist[["simul"]][[paste0(enh,"_conflow")]], enh_prop_dist[["simul"]][[paste0(enh,"_confup")]]))

ymax=ymax*1.1

par(mar=c(3.1, 4.5, 2.75, 1))

plot(enh_prop_dist[["obs"]][[enh]], col=dataset.colors["Original"], main="", type="n", xlab="",ylab="",  axes=F, ylim=c(ymin,ymax))

## lines(enh_prop_dist[["obs"]][[enh]], col=dataset.colors["Original"])
## lines(enh_prop_dist[["simul"]][[enh]], col=dataset.colors["Simulated"])

points(enh_prop_dist[["obs"]][[enh]], col=dataset.colors["Original"], pch=20)
points(enh_prop_dist[["simul"]][[enh]], col=dataset.colors["Simulated"], pch=20)

xpos=1:length(enh_prop_dist[["obs"]][[enh]])

segments(xpos, enh_prop_dist[["obs"]][[paste0(enh,"_conflow")]], xpos, enh_prop_dist[["obs"]][[paste0(enh,"_confup")]], col=dataset.colors["Original"])
segments(xpos, enh_prop_dist[["simul"]][[paste0(enh,"_conflow")]], xpos, enh_prop_dist[["simul"]][[paste0(enh,"_confup")]], col=dataset.colors["Simulated"])


class_leg <- c("0", "0.5", "1", "1.5", "2")
axis(side=1, at=c(1,10,20,30,40), labels=class_leg, mgp=c(3, 0.65, 0), cex.axis=1.1)
axis(side=2, mgp=c(3, 0.65, 0), cex.axis=1.1, las=2)

legend("topright", legend=c("PCHi-C data", "simulated data"), col=dataset.colors[c("Original", "Simulated")],lty=1, seg.len=1, bty='n', cex=1.1, inset=c(0.05, -0.1), xpd=NA)

mtext("% length covered by enhancers", side=2, cex=0.85, line=2.25, at=(ymin+ymax)*0.9/2)
mtext("distance to baits (Mb)", side=1, line=2.25, cex=0.85)

mtext("B", side=3, line=1.45, at=-5.75, font=2, cex=1.2)

#############################################################################################################################
##################################  Fig2-C - Enhancer proportion according to nb of cell types ##############################

ymax=max(c(enh_prop_nb_cell[["obs"]][[paste0(enh,"_conflow")]], enh_prop_nb_cell[["obs"]][[paste0(enh,"_confup")]],  enh_prop_nb_cell[["simul"]][[paste0(enh,"_conflow")]], enh_prop_nb_cell[["simul"]][[paste0(enh,"_confup")]]))

ylim=c(0, ymax)
xlim=c(0.5, length(enh_prop_nb_cell[["obs"]][[paste0(enh,"_conflow")]])+0.5)

par(mar=c(3.1, 4.5, 3, 1))

plot(enh_prop_nb_cell[["obs"]][[enh]], type="n", ylim=ylim,  xlim=xlim, axes=F, xlab="", ylab="", xaxs="i")

## lines(enh_prop_nb_cell[["obs"]][[enh]], col=dataset.colors["Original"])
## lines(enh_prop_nb_cell[["simul"]][[enh]], col=dataset.colors["Simulated"])

points(enh_prop_nb_cell[["obs"]][[enh]], col=dataset.colors["Original"], pch=20)
points(enh_prop_nb_cell[["simul"]][[enh]], col=dataset.colors["Simulated"], pch=20)

xpos=1:length(enh_prop_nb_cell[["obs"]][[enh]])

segments(xpos, enh_prop_nb_cell[["obs"]][[paste0(enh,"_conflow")]], xpos, enh_prop_nb_cell[["obs"]][[paste0(enh,"_confup")]], col=dataset.colors["Original"])
segments(xpos, enh_prop_nb_cell[["simul"]][[paste0(enh,"_conflow")]], xpos, enh_prop_nb_cell[["simul"]][[paste0(enh,"_confup")]], col=dataset.colors["Simulated"])

axis(side=1, mgp=c(3, 0.65, 0), cex.axis=1.1, at=seq(from=1, to=max(xpos), by=1))
axis(side=2, mgp=c(3, 0.65, 0), cex.axis=1.1, las=2)


mtext("% length covered by enhancers", side=2, cex=0.85, line=2.7, at=sum(ylim)*0.9/2)
mtext("number of cell types", side=1, line=2, cex=0.85)

at=-0.75
mtext("C", side=3, line=1.25, at=at, font=2, cex=1.2)

############################################################################################################################
############################################## correlation gene expression and enhancers activity ##########################

ymin=min(c(correl_activity[["obs"]][[paste0(enh,"_conflow")]], correl_activity[["obs"]][[paste0(enh,"_confup")]],  correl_activity[["simul"]][[paste0(enh,"_conflow")]], correl_activity[["simul"]][[paste0(enh,"_confup")]]))
ymax=max(c(correl_activity[["obs"]][[paste0(enh,"_conflow")]], correl_activity[["obs"]][[paste0(enh,"_confup")]],  correl_activity[["simul"]][[paste0(enh,"_conflow")]], correl_activity[["simul"]][[paste0(enh,"_confup")]]))

ylim=c(ymin, ymax)

par(mar=c(3.1, 4.5, 3, 1))

plot(correl_activity[["obs"]][[enh]], type="n", ylab="", main="", las=2, ylim=ylim, axes=F)

## lines(correl_activity[["obs"]][[enh]], col=dataset.colors["Original"])
## lines(correl_activity[["simul"]][[enh]], col=dataset.colors["Simulated"])

points(correl_activity[["obs"]][[enh]], col=dataset.colors["Original"], pch=20)
points(correl_activity[["simul"]][[enh]], col=dataset.colors["Simulated"], pch=20)


xpos=1:length(correl_activity[["obs"]][[enh]])

segments(xpos, correl_activity[["obs"]][[paste0(enh,"_conflow")]], xpos, correl_activity[["obs"]][[paste0(enh,"_confup")]], col=dataset.colors["Original"])
segments(xpos, correl_activity[["simul"]][[paste0(enh,"_conflow")]], xpos, correl_activity[["simul"]][[paste0(enh,"_confup")]], col=dataset.colors["Simulated"])



class_leg <- c("0", "0.5", "1", "1.5", "2")
axis(side=1, at=c(1,10,20,30,40), labels=class_leg, mgp=c(3, 0.65, 0), cex.axis=1.1)
axis(side=2, mgp=c(3, 0.65, 0), cex.axis=1.1, las=2)

mtext("Spearman's rho", side=2, cex=0.85, line=3, at=(ymin+ymax)*0.9/2)
mtext("distance to promoters (Mb)", side=1, line=2, cex=0.85)

mtext("D", side=3, line=1.25, at=-7.5, font=2, cex=1.2)

##################################################################################################################################

dev.off()


##################################################################################################################################
