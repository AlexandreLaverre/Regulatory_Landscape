#########################################################################################################################
## if it's the first time we run this figure, we load and prepare data

objects=ls()

if(!"pathScripts"%in%objects){
  load=T
  prepare=T
}

#########################################################################################################################

source("parameters.R") ## paths are defined based on the user name

#########################################################################################################################

if(load){
  ref_sp = "human"
  
  load(paste(pathFigures, "RData/Fig2_", ref_sp, "_A_B_C.Rdata", sep=""))
  load(paste(pathFigures, "RData/Fig2_", ref_sp, "_D_E.Rdata", sep=""))
  
  enhancers = enhancer.datasets[[ref_sp]]
}

#########################################################################################################################

pdf(paste(pathFigures, "Figure2.pdf", sep=""), width=8.5, height=5)

par(mai = c(0.5, 0.7, 0.3, 0.2)) # bottom, left, top, right
layout(matrix(c(1, 1, 2, 2, 3, 4, 5, 5), nrow = 2, byrow = TRUE))

############################################  Fig2-A - Global enhancer proportion ############################################ 

barcenter <- barplot(enh_prop$data, border=c(dataset.colors, "white"), col=c(dataset.colors, "white"),
                     lwd=1.5, cex.names=0.8, density=c(dataset.density,0), angle=c(dataset.angle,0),
                     ylim=c(0,15), ylab="Enhancer proportion (%)", axisnames = F, main="", las=2)


text(c(1.3,4.9,8.5,12.1), par("usr")[3]-0.005, labels = enh.syn, pos = 1, xpd = TRUE, cex=1)

legend("topleft", legend = c("Original", "Simulated"), 
       border=dataset.colors, density=c(dataset.density,0), angle=c(dataset.angle,0), bty='n', cex=1.2)
par(lwd=1)

segments(barcenter, enh_prop$conf_up, barcenter, enh_prop$conf_low, lwd = 3)
arrows(barcenter, enh_prop$conf_up, barcenter, enh_prop$conf_low, lwd = 1.5, angle = 90, code = 3, length = 0.05)

for (x in seq(1,length(barcenter)-1, by=3)){
  segments(barcenter[x], enh_prop$data[x]+1, barcenter[x+1], enh_prop$data[x]+1) 
  text("***", x=(barcenter[x]+barcenter[x+1])/2, y=enh_prop$data[x]+1.5)
}

mtext("A", side=3, line=1, at=-1.5, font=2, cex=1.2)

############################################  Fig2-B - Enhancer proportion according to distance ############################################ 
if(ref_sp=="human"){YMAX=0.2}else{YMAX=0.15}

plot(prop_dist[["obs"]]$FANTOM5, type="l", col="white", ylab="Enhancer length proportion (mean)", main="", las=2,
     xlab='', xaxt = "n", ylim=c(0,YMAX))

nb_col = 1
enhancers <- c("ENCODE", "FANTOM5")
if (ref_sp == "human"){enhancers <- c(enhancers, "RoadmapEpigenomics", "FOCS_GRO_seq")}

for (enh in enhancers){
  
  points(prop_dist[["obs"]][[paste0(enh)]], type="l", col=col.enhancers[nb_col])
  for (row in 1:nrow(prop_dist[["obs"]])){
    segments(x0=row,y0=prop_dist[["obs"]][row,paste0(enh, "_conflow")],
             x1=row,y1=prop_dist[["obs"]][row,paste0(enh, "_confup")], col=col.enhancers[nb_col], lwd=0.5)}
  
  # points(prop_dist[["simul"]][[paste0(enh)]][0:50], type="l", lty=2, col=col.enhancers[nb_col], lwd=0.6)
  # for (row in 1:nrow(prop_dist[["simul"]][0:50,])){
  #   segments(x0=row,y0=prop_dist[["simul"]][row,paste0(enh, "_conflow")],
  #            x1=row,y1=prop_dist[["simul"]][row,paste0(enh, "_confup")], col=col.enhancers[nb_col], lwd=0.5)}
  # 
  nb_col = nb_col + 1
  
}

class_leg <- c("0", "0.5", "1", "1.5", "2")
axis(side=1, at=c(1,10,20,30,40), labels=class_leg, mgp=c(3, 0.65, 0), cex.axis=1.1)
legend("topright", legend=label.enhancers, col=col.enhancers, bty='n', lty=1, cex=1, ncol=2)
mtext("B", side=3, line=1, at=-4.5, font=2, cex=1.2)
mtext("Linear distance to promoters regions (Mb)", side=1, line=2.25, cex=0.7)

############################################  Fig2-C - Enhancer proportion according to nb of sample ############################################ 
if(ref_sp=="human"){YMAX=0.2; x_leg_class=5}else{YMAX=0.15; x_leg_class=2}

plot(prop_nb_sample[["obs"]]$FANTOM5, col="white", ylim=c(0,YMAX), las=2,
     ylab="Enhancer length proportion (mean)", xlab="", xaxt = "n", main="")

col_nb = 1
for (enh in enhancers){
  points(prop_nb_sample[["obs"]][[enh]], type="l", col=col.enhancers[col_nb])
  for (row in 1:nrow(prop_nb_sample[["obs"]])){
    segments(x0=row,y0=prop_nb_sample[["obs"]][row,paste0(enh, "_conflow")],
             x1=row,y1=prop_nb_sample[["obs"]][row,paste0(enh, "_confup")], col=col.enhancers[col_nb], lwd=0.5)}
  
  # points(prop_nb_sample[["simul"]][[enh]], type="l", lty=3, col=col.enhancers[col_nb])
  # for (row in 1:nrow(prop_nb_sample[["simul"]])){
  #   segments(x0=row,y0=prop_nb_sample[["simul"]][row,paste0(enh, "_conflow")],
  #            x1=row,y1=prop_nb_sample[["simul"]][row,paste0(enh, "_confup")], col=col.enhancers[col_nb], lwd=0.3)}
  # 
  col_nb = col_nb + 1
}

if (ref_sp == "mouse"){class_leg <- c("1", "3", "5", "7", "9", "11", "13"); max_nb_cell=13
}else{class_leg <- c("1", "5", "10", "15", "20", "25"); max_nb_cell=25}

axis(side=1, at=seq(1,max_nb_cell+1, x_leg_class), labels=class_leg, mgp=c(3, 0.65, 0), cex.axis=1)

# axis(1, at=seq(1,max_nb_cell+1, x_leg_class), labels=F)
# text(seq(1,max_nb_cell+1, x_leg_class), par("usr")[3]-0.01, labels = class_leg, pos = 1, xpd = TRUE)
#legend("topleft", legend=enhancers, col=col.enhancers, bty='n', lty=1, cex=0.8)

mtext("C", side=3, line=1, at=-4.5, font=2, cex=1.2)
mtext("Number of sample", side=1, line=2.25, cex=0.7)


############################################ Fig2-D - Gene expression vs nb enhancers ############################################ 
plot(gene_expression_enhancers$ENCODE, type="l", col="white", ylab="Average expression level (log2 RPKM)", main="", las=2,
     xlab="", xaxt = "n", ylim=c(2.7,3.5))

nb_col = 1
for (enh in enhancers){
  
  points(gene_expression_enhancers[[paste0(enh)]], type="l", col=col.enhancers[nb_col])
  for (row in 1:nrow(gene_expression_enhancers)){
    segments(x0=row,y0=gene_expression_enhancers[row,paste0(enh, "_conflow")],
             x1=row,y1=gene_expression_enhancers[row,paste0(enh, "_confup")], col=col.enhancers[nb_col], lwd=0.5)}

  nb_col = nb_col + 1
  
}


axis(side=1, at=seq(1,10,1), labels=seq(1,10,1), mgp=c(3, 0.65, 0), cex.axis=1)

mtext("D", side=3, line=1, at=-1, font=2, cex=1.2)
mtext("Quantile of Number of contacted enhancers", side=1, line=2.25, cex=0.7)

############################################ Fig2-E - Correlation Gene expression and enhancers activity ############################################ 
if(ref_sp=="human"){YMAX=0.2}else{YMAX=0.15}

plot(correl_activity[["obs"]]$FANTOM5, type="l", col="white", ylab="Spearman's correlation coefficient (mean)", main="", las=2,
     xlab="", xaxt = "n", ylim=c(0,0.35))

nb_col = 1
for (enh in enhancers){
  
  points(correl_activity[["obs"]][[paste0(enh)]], type="l", col=col.enhancers[nb_col])
  for (row in 1:nrow(correl_activity[["obs"]])){
    segments(x0=row,y0=correl_activity[["obs"]][row,paste0(enh, "_conflow")],
             x1=row,y1=correl_activity[["obs"]][row,paste0(enh, "_confup")], col=col.enhancers[nb_col], lwd=0.5)}
  
  points(correl_activity[["simul"]][[paste0(enh)]], type="l", lty=2, col=col.enhancers[nb_col], lwd=0.6)
  for (row in 1:nrow(correl_activity[["simul"]])){
    segments(x0=row,y0=correl_activity[["simul"]][row,paste0(enh, "_conflow")],
             x1=row,y1=correl_activity[["simul"]][row,paste0(enh, "_confup")], col=col.enhancers[nb_col], lwd=0.5)}

  nb_col = nb_col + 1
  
}

class_leg <- c("0", "0.5", "1", "1.5", "2")
axis(side=1, at=c(1,10,20,30,40), labels=class_leg, mgp=c(3, 0.65, 0), cex.axis=1.1)
legend("topright", legend="Simulated", col="black", bty='n', lty=2, cex=1.2)
mtext("Linear distance to promoters regions (Mb)", side=1, line=2.25, cex=0.7)

mtext("E", side=3, line=1, at=-4.5, font=2, cex=1.2)

dev.off()

