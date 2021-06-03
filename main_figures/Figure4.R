#########################################################################################################################

objects=ls()

if(!"pathFigures"%in%objects){

  source("parameters.R")
  
  ref_sp <- "human"
  target_sp = setdiff(c("human", "mouse"), ref_sp)

  if(ref_sp == "human"){
    species <-c("macaque", "mouse", "rat",  "rabbit", "dog", "cow", "elephant", "opossum", "chicken")
  }
  
  if(ref_sp == "mouse"){
    species <-c("rat", "rabbit", "macaque", "human", "dog", "cow", "elephant", "opossum", "chicken")
  } 
    
  load=T
  prepare=T
}

#########################################################################################################################

if(load){
  print("loading data")
  
  load(paste(pathFigures, "RData/data.synteny.conservation.",ref_sp,".RData", sep=""))
  
  load=FALSE
}

#########################################################################################################################

if(prepare){
  enh="ENCODE"

  ## compute statistics for all species, all interactions between minDistanceSyntenyRef and maxDistanceSyntenyRef
  
  prop.obs.alldist <- c()
  prop.sim.alldist <- c()

  ci.low.obs.alldist <- c()
  ci.low.sim.alldist <- c()
  
  ci.high.obs.alldist <- c()
  ci.high.sim.alldist <- c()

  pvalues.alldist <- c()
    
  for(sp in species){
    synt_obs <- conserv_synteny[[enh]][[sp]][["synt_obs"]]
    synt_simul <- conserv_synteny[[enh]][[sp]][["synt_simul"]]
    
    nb_cons_synt_obs <- length(which(synt_obs$origin_dist >= minDistanceSyntenyRef & synt_obs$origin_dist <= maxDistanceSyntenyRef & synt_obs$target_dist <= maxDistanceSyntenyTarget))
    nb_cons_synt_simul <- length(which(synt_simul$origin_dist >= minDistanceSyntenyRef & synt_simul$origin_dist <= maxDistanceSyntenyRef & synt_simul$target_dist <= maxDistanceSyntenyTarget))
    
    nb_tot_synt_obs <- length(which(synt_obs$origin_dist >= minDistanceSyntenyRef & synt_obs$origin_dist <= maxDistanceSyntenyRef))
    nb_tot_synt_simul <- length(which(synt_simul$origin_dist >= minDistanceSyntenyRef & synt_simul$origin_dist <= maxDistanceSyntenyRef))

    ## chi-squared test
    mat <- matrix(c(nb_cons_synt_obs, nb_cons_synt_simul, nb_tot_synt_obs-nb_cons_synt_obs, nb_tot_synt_simul-nb_cons_synt_simul), nrow=2)
    pval <-chisq.test(mat)$p.value
    
    pvalues.alldist <- c(pvalues.alldist, pval)

    prop.obs.alldist <- c(prop.obs.alldist, 100*nb_cons_synt_obs/nb_tot_synt_obs)
    prop.sim.alldist <- c(prop.sim.alldist, 100*nb_cons_synt_simul/nb_tot_synt_simul)
    
    prop.test.obs <- prop.test(x = nb_cons_synt_obs, n=nb_tot_synt_obs, p=0.5)
    prop.test.simul <- prop.test(x = nb_cons_synt_simul, n=nb_tot_synt_simul, p=0.5)

    ci.low.obs.alldist <- c(ci.low.obs.alldist, 100*prop.test.obs$conf.int[1])
    ci.high.obs.alldist <- c(ci.high.obs.alldist, 100*prop.test.obs$conf.int[2])

    ci.low.sim.alldist <- c(ci.low.sim.alldist, 100*prop.test.simul$conf.int[1])
    ci.high.sim.alldist <- c(ci.high.sim.alldist, 100*prop.test.simul$conf.int[2])
    
    print(paste0(ref_sp, " to ", sp))
    print(paste0("PCHi-C promoters-ENCODE pairs maintained in synteny : ", round(100*nb_cons_synt_obs/nb_tot_synt_obs,2), "%"))
    print(paste0("Simulated promoters-ENCODE pairs maintained in synteny : ", round(100*nb_cons_synt_simul/nb_tot_synt_simul,2), "%"))
    print(paste0("Wilcoxon test p-value: ", pval))
  }
  
  names(prop.obs.alldist) <- species
  names(prop.sim.alldist) <- species

  names(ci.low.obs.alldist) <- species
  names(ci.high.obs.alldist) <- species
  
  names(ci.low.sim.alldist) <- species
  names(ci.high.sim.alldist) <- species

  names(pvalues.alldist) <- species

  
  prepare=FALSE
}

#########################################################################################################################

## 1 column width 85 mm = 3.34 in
## 1.5 column width 114 mm = 4.49 in
## 2 columns width 174 mm = 6.85 in
## max height: 11 in

#########################################################################################################################

pdf(paste(pathFigures, "Figure4.pdf", sep=""), width=6.85, height=5.5)

layout(matrix(c(1, 1, 2, 3), nrow = 2, byrow = TRUE))

#########################################################################################################################

## global synteny conservation, ENCODE

enh="ENCODE"

m <- matrix(c(prop.obs.alldist, prop.sim.alldist), nrow=2, byrow=T)

par(mar=c(3.1, 4.1, 2.1, 1.1))

bar <-barplot(m, beside=T, space=c(0.25, 1.2), col=dataset.colors, border=dataset.colors, axes=F,  ylim = c(55,105), xpd=F)
colnames(bar)<-species

axis(side=2, las=2,  mgp=c(3, 0.75, 0), cex.axis=0.9)
mtext("% pairs in conserved synteny", side=2, line=2.5, cex=0.8)

xax=apply(bar, 2, mean)
axis(side=1, at=xax, mgp=c(3, 0.65, 0), labels=rep("", length(species)))
mtext(species, side=1, at=xax, line=0.5, cex=0.75)

mtext(paste(ref_sp, "vs."), side=1, line=0.5, at=min(xax)-diff(xax)[1]*1, cex=0.75)

segments(bar[1,], ci.low.obs.alldist, bar[1,], ci.high.obs.alldist, col="black", lwd=1.1)

segments(bar[2,], ci.low.sim.alldist, bar[2,], ci.high.sim.alldist, col="black", lwd=1.1)


smallx=(bar[2,2]-bar[1,2])/10

for (sp in species){
  this.pval=pvalues.alldist[sp]

  if(this.pval < 0.0001){
    text="***"
  } else{
    if(this.pval < 0.001){
      text="**"
    } else{
      if(this.pval < 0.01){
        text="*"
      } else{
        text="NS"
      }
    }
  }

  ypos=max(prop.obs.alldist[sp], prop.sim.alldist[sp])+2
  
  segments(bar[1,sp]+smallx, ypos, bar[2,sp]-smallx, ypos)
  text(text, x=mean(as.numeric(bar[,sp])), y=ypos+2.5, xpd=NA, cex=1)
}

## legend
legend("topright", legend = c("PCHi-C data", "simulated data"), fill=dataset.colors, border=dataset.colors, bty='n', cex=0.95, inset=c(0, -0.15), xpd=NA)

## plot label
mtext("A", side=3, line=1, at=xax[1]-diff(xax)[1]*1.42, font=2, cex=1.005)

###################################################################################

## synteny conservation as a function of the distance between promoter and enhancer

selected.species <- c("mouse", "opossum")
labels <- c("B", "C")
names(labels) <- selected.species

nbkept=25

for (sp in selected.species){
  conserv.obs <- conserv_synteny[[enh]][[sp]][["synt_obs"]]
  conserv.sim <- conserv_synteny[[enh]][[sp]][["synt_simul"]]

  prop.obs=tapply(conserv.obs$target_dist, conserv.obs$class_dist, function(x) 100*length(which(x<maxDistanceSyntenyTarget))/length(x))[1:nbkept]
  prop.sim=tapply(conserv.sim$target_dist, conserv.sim$class_dist, function(x) 100*length(which(x<maxDistanceSyntenyTarget))/length(x))[1:nbkept]

  ci.obs.low=tapply(conserv.obs$target_dist, conserv.obs$class_dist, function(x) 100*prop.test(length(which(x<maxDistanceSyntenyTarget)), length(x))$conf.int[1])[1:nbkept]
  ci.obs.high=tapply(conserv.obs$target_dist, conserv.obs$class_dist, function(x) 100*prop.test(length(which(x<maxDistanceSyntenyTarget)), length(x))$conf.int[2])[1:nbkept]

  ci.sim.low=tapply(conserv.sim$target_dist, conserv.sim$class_dist, function(x) 100*prop.test(length(which(x<maxDistanceSyntenyTarget)), length(x))$conf.int[1])[1:nbkept]
  ci.sim.high=tapply(conserv.sim$target_dist, conserv.sim$class_dist, function(x) 100*prop.test(length(which(x<maxDistanceSyntenyTarget)), length(x))$conf.int[2])[1:nbkept]
  
  par(mar=c(3.1, 4.1, 2.1, 1.1))
  
  ylim=range(c(ci.sim.low, ci.sim.high, ci.obs.low, ci.obs.high))

  smally=diff(ylim)/10
  ylim=ylim+c(-smally, smally)
  
  
  plot(1, type="n", xlim=c(0.5, nbkept+0.5), ylim=ylim, xlab="", ylab="", axes=F)

  xpos=1:nbkept

  ## lines(xpos, prop.obs, col=dataset.colors["Original"])
  points(xpos, prop.obs, col=dataset.colors["Original"], pch=20)
  segments(xpos, ci.obs.low, xpos, ci.obs.high, col=dataset.colors["Original"])

  ## lines(xpos, prop.sim, col=dataset.colors["Simulated"])
  points(xpos, prop.sim, col=dataset.colors["Simulated"], pch=20)
  segments(xpos, ci.sim.low, xpos, ci.sim.high, col=dataset.colors["Simulated"])

  axis(side=2, las=2,  mgp=c(3, 0.75, 0), cex.axis=0.9)
  mtext("% pairs in conserved synteny", side=2, line=2.5, cex=0.8)
  
  xax=c(1, 5, 10, 15, 20, 25)
  axlab=as.character(c(0.05, 0.25, 0.5, 0.75, 1, 1.5))
  
  axis(side=1, mgp=c(3, 0.65, 0), at=xax, labels=axlab, cex.axis=0.9)
  mtext("distance to promoters (Mb)", side=1, line=2, cex=0.8)

  mtext(paste(ref_sp, "vs.", sp, sep=" "), side=3, cex=0.8, line=0.5)


  ## plot label
  mtext(labels[sp], side=3, line=1, at=-5.5, font=2, cex=1.005)
  
}

#####################################################################################

dev.off()

#######################################################################################
