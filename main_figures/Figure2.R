#################################################################################################################
## if it's the first time we run this figure, we load and prepare data

objects=ls()

if(!"pathScripts"%in%objects){
  load=T
  prepare=T
}

source("parameters.R") ## paths are defined based on the user name

#################################################################################################################

if(load){
  ref_sp = "human"
  
  load(paste(pathFigures, "RData/data.enhancer.coverage.", ref_sp, ".Rdata", sep=""))
  load(paste(pathFigures, "RData/data.promoter.enhancer.correlation.", ref_sp, ".Rdata", sep=""))
  
  enhancers = enhancer.datasets[[ref_sp]]
}

#################################################################################################################

## 1 column width 85 mm = 3.34 in
## 1.5 column width 114 mm = 4.49 in 
## 2 columns width 174 mm = 6.85 in
## max height: 11 in

#################################################################################################################

pdf(paste(pathFigures, "Figure2.pdf", sep=""), width=6.85, height=5)

par(mai = c(0.5, 0.5, 0.3, 0.2)) # bottom, left, top, right

layout(matrix(c(1, 1, 2, 2, 3, 4, 5, 5), nrow = 2, byrow = TRUE))

##################################################################################################################
############################################  % length covered by enhancers  #####################################

par(mar=c(4.1, 4.5, 2.75, 1))

m.prop=t(matrix(enh_prop$data, nrow=4, byrow=T))

barcenter <- barplot(m.prop, beside=T,  border=dataset.colors, col=dataset.colors, 
                     lwd=1.5, cex.names=0.8, density=dataset.density, angle=dataset.angle,
                     ylim=c(0,15), ylab="", xlab="", axisnames = F, main="", space=c(0.1, 1), las=2)


xposlab=apply(barcenter, 2, mean)
allxpos=as.numeric(barcenter)

mtext(enh.syn.narrow, line=c(rep(0.5,2), rep(1.3,2)), side=1, at=xposlab, cex=0.75)

mtext("% length covered by enhancers", side=2, cex=0.85, line=2.7, at=7)

legend("topleft", legend = c("PCHi-C data", "simulated data"), border=dataset.colors, density=c(dataset.density), angle=c(dataset.angle), bty='n', cex=1.1, inset=c(0, -0.1), xpd=NA)

par(lwd=1)
segments(allxpos, enh_prop$conf_up, allxpos, enh_prop$conf_low, lwd = 3)
arrows(allxpos, enh_prop$conf_up, allxpos, enh_prop$conf_low, lwd = 1.5, angle = 90, code = 3, length = 0.05)

for (x in seq(1,length(allxpos)-1, by=2)){
   segments(allxpos[x], enh_prop$data[x]+1, allxpos[x+1], enh_prop$data[x]+1) 
   text("***", x=(allxpos[x]+allxpos[x+1])/2, y=enh_prop$data[x]+1.5, cex=1.2)
 }

mtext("a", side=3, line=2, at=-1.4, font=2, cex=1.2)

###########################################################################################################################
##########################  Fig2-B - Enhancer proportion according to distance ############################################

## only ENCODE

enh="ENCODE"

ymax=max(c(enh_prop_dist[["obs"]][[paste0(enh,"_conflow")]], enh_prop_dist[["obs"]][[paste0(enh,"_confup")]],  enh_prop_dist[["simul"]][[paste0(enh,"_conflow")]], enh_prop_dist[["simul"]][[paste0(enh,"_confup")]]))
ymin=min(c(enh_prop_dist[["obs"]][[paste0(enh,"_conflow")]], enh_prop_dist[["obs"]][[paste0(enh,"_confup")]],  enh_prop_dist[["simul"]][[paste0(enh,"_conflow")]], enh_prop_dist[["simul"]][[paste0(enh,"_confup")]]))

ymax=ymax*1.1

plot(enh_prop_dist[["obs"]][[enh]], col=dataset.colors["Original"], main="", type="n", xlab="",ylab="",  axes=F, ylim=c(ymin,ymax))

lines(enh_prop_dist[["obs"]][[enh]], col=dataset.colors["Original"])
lines(enh_prop_dist[["simul"]][[enh]], col=dataset.colors["Simulated"])


xpos=1:length(enh_prop_dist[["obs"]][[enh]])

segments(xpos, enh_prop_dist[["obs"]][[paste0(enh,"_conflow")]], xpos, enh_prop_dist[["obs"]][[paste0(enh,"_confup")]], col=dataset.colors["Original"])
segments(xpos, enh_prop_dist[["simul"]][[paste0(enh,"_conflow")]], xpos, enh_prop_dist[["simul"]][[paste0(enh,"_confup")]], col=dataset.colors["Simuated"])


class_leg <- c("0", "0.5", "1", "1.5", "2")
axis(side=1, at=c(1,10,20,30,40), labels=class_leg, mgp=c(3, 0.65, 0), cex.axis=1.1)
axis(side=2, mgp=c(3, 0.65, 0), cex.axis=1.1, las=2)

legend("topright", legend=c("PCHi-C data", "simulated data"), col=dataset.colors[c("Original", "Simulated")],lty=1, seg.len=1, bty='n', cex=1.1, inset=c(0.05, -0.1), xpd=NA)

mtext("% length covered by enhancers", side=2, cex=0.85, line=2.7, at=(ymin+ymax)*0.8/2)
mtext("linear distance to promoter regions (Mb)", side=1, line=2.25, cex=0.85)

mtext("b", side=3, line=1, at=-5, font=2, cex=1.2)

## ###########################################################################################################################
## ################################  Fig2-C - Enhancer proportion according to nb of cell types ##############################
## if(ref_sp=="human"){YMAX=0.2; x_leg_class=5}else{YMAX=0.15; x_leg_class=2}

## plot(prop_nb_sample[["obs"]]$FANTOM5, col="white", ylim=c(0,YMAX), las=2,
##      ylab="Enhancer length proportion (mean)", xlab="", xaxt = "n", main="")

## col_nb = 1
## for (enh in enhancers){
##   points(prop_nb_sample[["obs"]][[enh]], type="l", col=col.enhancers[col_nb])
##   for (row in 1:nrow(prop_nb_sample[["obs"]])){
##     segments(x0=row,y0=prop_nb_sample[["obs"]][row,paste0(enh, "_conflow")],
##              x1=row,y1=prop_nb_sample[["obs"]][row,paste0(enh, "_confup")], col=col.enhancers[col_nb], lwd=0.5)}
  
##   col_nb = col_nb + 1
## }

## if (ref_sp == "mouse"){class_leg <- c("1", "3", "5", "7", "9", "11", "13"); max_nb_cell=13
## }else{class_leg <- c("1", "5", "10", "15", "20", "25"); max_nb_cell=25}

## at=seq(0,max_nb_cell+1, x_leg_class)
## at[1] <- 1
## axis(side=1, at=at, labels=class_leg, mgp=c(3, 0.65, 0), cex.axis=1)

## mtext("C", side=3, line=1, at=-4.5, font=2, cex=1.2)
## mtext("Number of cell types", side=1, line=2.25, cex=0.7)

## ##################################################################################################################################
## ############################################ Fig2-D - Gene expression vs nb enhancers ############################################ 
## plot(gene_expression_enhancers$ENCODE, type="l", col="white", ylab="Average expression level (log2 RPKM)", main="", las=2,
##      xlab="", xaxt = "n", ylim=c(2.7,3.5))

## nb_col = 1
## for (enh in enhancers){
  
##   points(gene_expression_enhancers[[paste0(enh)]], type="l", col=col.enhancers[nb_col])
##   for (row in 1:nrow(gene_expression_enhancers)){
##     segments(x0=row,y0=gene_expression_enhancers[row,paste0(enh, "_conflow")],
##              x1=row,y1=gene_expression_enhancers[row,paste0(enh, "_confup")], col=col.enhancers[nb_col], lwd=0.5)}

##   nb_col = nb_col + 1
## }

## axis(side=1, at=seq(1,10,1), labels=seq(1,10,1), mgp=c(3, 0.65, 0), cex.axis=1)

## mtext("D", side=3, line=1, at=-1, font=2, cex=1.2)
## mtext("Quantile of Number of contacted enhancers", side=1, line=2.25, cex=0.7)

## ##################################################################################################################################
## ############################################ Fig2-E - Correlation Gene expression and enhancers activity ##########################
## if(ref_sp=="human"){YMAX=0.2}else{YMAX=0.15}

## plot(correl_activity[["obs"]]$FANTOM5, type="l", col="white", ylab="Spearman's correlation coefficient (mean)", main="", las=2,
##      xlab="", xaxt = "n", ylim=c(0,0.35))

## nb_col = 1
## for (enh in enhancers){
  
##   points(correl_activity[["obs"]][[paste0(enh)]], type="l", col=col.enhancers[nb_col])
##   for (row in 1:nrow(correl_activity[["obs"]])){
##     segments(x0=row,y0=correl_activity[["obs"]][row,paste0(enh, "_conflow")],
##              x1=row,y1=correl_activity[["obs"]][row,paste0(enh, "_confup")], col=col.enhancers[nb_col], lwd=0.5)}
  
##   points(correl_activity[["simul"]][[paste0(enh)]], type="l", lty=2, col=col.enhancers[nb_col], lwd=0.6)
##   for (row in 1:nrow(correl_activity[["simul"]])){
##     segments(x0=row,y0=correl_activity[["simul"]][row,paste0(enh, "_conflow")],
##              x1=row,y1=correl_activity[["simul"]][row,paste0(enh, "_confup")], col=col.enhancers[nb_col], lwd=0.5)}

##   nb_col = nb_col + 1
## }

## class_leg <- c("0", "0.5", "1", "1.5", "2")
## axis(side=1, at=c(1,10,20,30,40), labels=class_leg, mgp=c(3, 0.65, 0), cex.axis=1.1)
## legend("topright", legend="Simulated", col="black", bty='n', lty=2, cex=1.2)
## mtext("Linear distance to promoters regions (Mb)", side=1, line=2.25, cex=0.7)

## mtext("E", side=3, line=1, at=-4.5, font=2, cex=1.2)

dev.off()


##################################################################################################################################
