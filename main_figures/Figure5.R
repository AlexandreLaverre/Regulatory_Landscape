##############################################################################

objects=ls()

if(!"pathScripts"%in%objects){
  load=T
  prepare=T

  source("parameters.R")
}

##############################################################################

if(load){
  ref_sp = "human"
  target_sp = "mouse"
  
  enhancers=enhancer.datasets[[ref_sp]]

  load(paste(pathFigures, "RData/data.sample.info.RData", sep=""))
  
  ## contact conservation data, already filtered to keep only filtered enhancers
  
  load(paste(pathFigures, "RData/data.contact.conservation.enhancers.RData", sep=""))

  contact.cons=contact.conservation[[paste0(ref_sp, "2", target_sp)]]

  sampleinfo.ref=sampleinfo[[ref_sp]]
  sampleinfo.tg=sampleinfo[[target_sp]]

  
  load=FALSE
}

##############################################################################

if(prepare){
  ## % conservation and confidence intervals, all enhancer datasets

  cons.stats <- list()

  for(enh in enhancers){

    ## contact conservation
    cc.obs <- contact.cons[[enh]][["obsobs"]]
    cc.sim <- contact.cons[[enh]][["simsim"]]
    
    cc.obs$dist_class <- cut(cc.obs$origin_dist, breaks=seq(from=minDistance, to=maxDistance, by=50e3), include.lowest = T)

    cc.obs$cons=apply(cc.obs[,sampleinfo.tg$Sample.ID], 1, function(x) any(x>0))
    cc.sim$cons=apply(cc.sim[,sampleinfo.tg$Sample.ID], 1, function(x) any(x>0))
    
    pc.cons.obs <- 100*length(which(cc.obs$cons))/dim(cc.obs)[1]
    pc.cons.sim <- 100*length(which(cc.sim$cons))/dim(cc.sim)[1]
    
    ## divide interactions based on distance, compute mean number of samples by distance class
    log((conserv$result-conserv$simul)/conserv$simul)
    
    pc.cons.dist <- sapply(as.factor(cc.obs$dist_class), function(x) 100*length(which(cc.obs$cons))/dim(cc.obs)[1]
                                                                 
    mean_dist <- t(sapply(filtered_data, function(x) tapply(x$distance, as.factor(x$dist_class), mean, na.rm=T)))
    
    dist_conf_low_celltypes <- t(sapply(filtered_data, function(x) tapply(x$nb_celltypes, as.factor(x$dist_class), function(y) {z<-t.test(y); return(z[["conf.int"]][1])})))
    dist_conf_high_celltypes <- t(sapply(filtered_data, function(x) tapply(x$nb_celltypes, as.factor(x$dist_class), function(y) {z<-t.test(y); return(z[["conf.int"]][2])})))
    

    test.obs <- prop.test(length(which(cons.obs)), dim(cc.obs)[1])
    test.sim <- prop.test(length(which(cons.sim)), dim(cc.sim)[1])
    
    cons.stats[[enh]]=list("cons.obs"=cons.obs, "cons.sim"=cons.sim,
                           "pc.cons.obs"=pc.cons.obs, "pc.cons.sim"=pc.cons.sim,
                           "test.obs"=test.obs, "test.sim"=test.sim)
    
  }
  
  prepare=FALSE
}

##############################################################################

## 1 column width 85 mm = 3.34 in
## 1.5 column width 114 mm = 4.49 in
## 2 columns width 174 mm = 6.85 in
## max height: 11 in

##############################################################################

#pdf(paste(pathFigures, "Figure5.pdf", sep=""), width=8.5, height=8)

m=matrix(rep(NA, 2*10), nrow=2)
m[1,]=c(rep(1,3), rep(2,2), rep(3,5))
m[2,]=c(rep(4,5), rep(5,5))

layout(m)

#################### Fig 5.A - % of conserved contacts #####################

if (ref_sp == "human"){YMAX=40}else{YMAX=30}

par(lwd = 1.5)

cons = sapply(enhancers, function(x) c(cons.stats[[x]]$pc.cons.obs, cons.stats[[x]]$pc.cons.sim))
conf.low =  sapply(enhancers, function(x) c(cons.stats[[x]]$test.obs$conf.int[1]*100, cons.stats[[x]]$test.sim$conf.int[1]*100))
conf.high =  sapply(enhancers, function(x) c(cons.stats[[x]]$test.obs$conf.int[2]*100, cons.stats[[x]]$test.sim$conf.int[2]*100))

b=barplot(cons, beside=T, names=rep("", dim(cons)[2]), ylim=c(0,YMAX), space=c(0.2,1),
          border=dataset.colors[c("Original", "Simulated")],  col=dataset.colors[c("Original", "Simulated")],
          mgp=c(3, 0.75, 0), cex.axis=1.1, las=2)


arrows(x0=b,y0=conf.low,y1=conf.high,angle=90,code=3,length=0.05)

## axis labels
mtext("% of conserved contacts", side=2, line=2.5, cex=1)
mtext(label.enhancers, at=apply(b, 2, mean), side=1, line=0.5, cex=1)

## legend & plot label
legend("topright", legend=c("PCHi-C data", "simulated data"), border=dataset.colors[c("Original", "Simulated")],
       fill=dataset.colors[c("Original", "Simulated")], bty='n', 
       cex=1.1, inset=c(0.05, -0.1), xpd=NA)

mtext("a", side=3, line=1, at=0.1, font=2, cex=1.2)

## mtext("A", side=3, line=1, at=0.1, font=2, cex=1.2)

## ############################################  B - Contact conservation by distance from TSS ############################################
## if (ref_sp == "human"){YLIM=c(-2,7)}else{YLIM=c(-5, 20)}
## class_leg <- c("0",  "0.5",  "1", "1.5", "2")

## par(lwd = 0.7)
## for (enh in enhancer.datasets[[ref_sp]]){
##   conserv = get(paste("conserv_dist", enh, sep="_"))
##   if (enh == "ENCODE"){ # First enhancers dataset
##     plot(log((conserv$result-conserv$simul)/conserv$simul), pch=20, col=col.enhancers[[enh]], xaxt = "n", ylim=YLIM,
##          xlab="Distance from TSS (Mb)", ylab="Excess of contact conservation\n from expected (log)", main="", las=2)
##   }else{points(log((conserv$result-conserv$simul)/conserv$simul), pch=20, col=col.enhancers[enh])} # Add lines of other enhancers datasets

##   for (row in 1:nrow(conserv)){
##     segments(x0=row,y0=log((conserv[row,]$obs_conflow-conserv[row,]$simul_conflow)/conserv[row,]$simul_conflow),
##              x1=row,y1=log((conserv[row,]$obs_confup-conserv[row,]$simul_confup)/conserv[row,]$simul_confup),
##              col=col.enhancers[enh], lwd=0.3)}

## }


## axis(side=1, at=c(1,10,20,30,40), labels=class_leg, mgp=c(3, 0.65, 0), cex.axis=1.1)
## legend("topleft", col=col.enhancers, legend = label.enhancers, bty='n',pch=20, cex=1)
## mtext("B", side=3, line=1, at=0.1, font=2, cex=1.2)

## ############################################  C - Contact conservation by nb samples ############################################
## ylim=c(0, 40)
## xlim=c(0.5, 8.5)
## xpos=seq(1, 8, 1)
## names(xpos) = 1:8

## smallx=c(-0.15, -0.075, 0.075, 0.15)
## names(smallx)=enhancer.datasets[[ref_sp]]

## plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=ylim, xaxs="i", yaxs="i")

## for(enh in enhancer.datasets[[ref_sp]]){
##   conserv = get(paste("conserv_sample", enh, sep="_"))
##   conserv = conserv[which(conserv$data == "obs"),]

##   for(row in 1:nrow(conserv)){

##     x=xpos[row]+smallx[enh]
##     points(x, conserv[row,'result'], pch=20, col=col.enhancers[enh], cex=1.1)
##     segments(x, conserv[row,'conflow'], x, conserv[row,'confup'], col=col.enhancers[enh])
##   }
## }

## abline(v=xpos[1:7]+0.5, lty=3, col="gray40")
## axis(side=1, at=xpos, mgp=c(3, 0.5, 0), labels=rep("", 8))
## mtext(conserv$class, at=xpos, side=1, line=1)
## mtext("Number of cell types", side=1, line=2.5)

## axis(side=2, mgp=c(3, 0.75, 0))
## mtext("Conserved contact (%)", side=2, line=2.5)

## legend("topleft", legend=label.enhancers, pch=20,
##        col=col.enhancers, bty="n", cex=1)

## mtext("C", side=3, line=1, at=0.1, font=2, cex=1.2)


## ############################################ D - Contact conservation by similar samples ####################################
## if (ref_sp == "human"){YMAX=30; side="topright"}else{YMAX=17; side="topleft"}

## par(lwd = 1.5)
## enh="ENCODE"
## conserv = get(paste("conserv_similar_sample", enh, sep="_"))

## bar <- barplot(conserv$data, border=c(dataset.colors,"white"), col=c(dataset.colors,"white"),
##                cex.names=0.8, density=c(dataset.density,0), angle=c(dataset.angle,0),
##                beside=T, space = c(0, 0.1, 0), las=2,
##                ylim=c(0,YMAX), main="", ylab="Conserved contact (%)")

## arrows(x0=bar,y0=conserv$conf_up,y1=conserv$conf_low,angle=90,code=3,length=0.05)
## #text(conserv$n_total, x=bar, y=conserv$conf_up+2, cex=0.8)

## cells = c("Pre-adipocytes", "ESC", "Bcell")
## text(c(1,4.2,7.4), par("usr")[3]-0.5, labels = cells, pos = 1, xpd = TRUE)
## legend(side, legend = c("Original", "Simulated"),
##        border=dataset.colors, fill=dataset.colors, density=dataset.density, angle=dataset.angle, bty='n')
## mtext("D", side=3, line=1, at=0.1, font=2, cex=1.2)

####################################################################################

#dev.off()

####################################################################################

