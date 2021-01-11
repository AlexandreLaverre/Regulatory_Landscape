#########################################################################################################################

objects=ls()

if(!"pathFigures"%in%objects){
  
  source("../main_figures/parameters.R")
  
  ref_sp <- "mouse"
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
  prop.obs.alldist <- list()
  prop.sim.alldist <- list()
  ci.low.obs.alldist <- list()
  ci.low.sim.alldist <- list()
  ci.high.obs.alldist <- list()
  ci.high.sim.alldist <- list()
  pvalues.alldist <- list()
  
  ## compute statistics for all enhancers, all species, all interactions between minDistanceSyntenyRef and maxDistanceSyntenyRef
  for (enh in enhancer.datasets[[ref_sp]]){
    
    for(sp in species){
      synt_obs <- conserv_synteny[[enh]][[sp]][["synt_obs"]]
      synt_simul <- conserv_synteny[[enh]][[sp]][["synt_simul"]]
      
      nb_cons_synt_obs <- length(which(synt_obs$origin_dist >= minDistanceSyntenyRef & synt_obs$origin_dist <= maxDistanceSyntenyRef & synt_obs$target_dist <= maxDistanceSyntenyTarget))
      nb_cons_synt_simul <- length(which(synt_simul$origin_dist >= minDistanceSyntenyRef & synt_simul$origin_dist <= maxDistanceSyntenyRef & synt_simul$target_dist <= maxDistanceSyntenyTarget))
      
      nb_tot_synt_obs <- length(which(synt_obs$origin_dist >= minDistanceSyntenyRef & synt_obs$origin_dist <= maxDistanceSyntenyRef))
      nb_tot_synt_simul <- length(which(synt_simul$origin_dist >= minDistanceSyntenyRef & synt_simul$origin_dist <= maxDistanceSyntenyRef))
      
      ## chi-squared test
      mat <- matrix(c(nb_cons_synt_obs, nb_cons_synt_simul, nb_tot_synt_obs-nb_cons_synt_obs, nb_tot_synt_simul-nb_cons_synt_simul), nrow=2)
      pval <- chisq.test(mat)$p.value
      
      pvalues.alldist[[enh]] <- c(pvalues.alldist[[enh]], pval)
      
      prop.obs.alldist[[enh]] <- c(prop.obs.alldist[[enh]], 100*nb_cons_synt_obs/nb_tot_synt_obs)
      prop.sim.alldist[[enh]] <- c(prop.sim.alldist[[enh]], 100*nb_cons_synt_simul/nb_tot_synt_simul)
      
      prop.test.obs <- prop.test(x = nb_cons_synt_obs, n=nb_tot_synt_obs, p=0.5)
      prop.test.simul <- prop.test(x = nb_cons_synt_simul, n=nb_tot_synt_simul, p=0.5)
      
      ci.low.obs.alldist[[enh]] <- c(ci.low.obs.alldist[[enh]], 100*prop.test.obs$conf.int[1])
      ci.high.obs.alldist[[enh]] <- c(ci.high.obs.alldist[[enh]], 100*prop.test.obs$conf.int[2])
      
      ci.low.sim.alldist[[enh]] <- c(ci.low.sim.alldist[[enh]], 100*prop.test.simul$conf.int[1])
      ci.high.sim.alldist[[enh]] <- c(ci.high.sim.alldist[[enh]], 100*prop.test.simul$conf.int[2])
    }
    
    names(prop.obs.alldist[[enh]]) <- species
    names(prop.sim.alldist[[enh]]) <- species
    
    names(ci.low.obs.alldist[[enh]]) <- species
    names(ci.high.obs.alldist[[enh]]) <- species
    
    names(ci.low.sim.alldist[[enh]]) <- species
    names(ci.high.sim.alldist[[enh]]) <- species
    
    names(pvalues.alldist[[enh]]) <- species
    
  }
  
  prepare=FALSE
}

#########################################################################################################################

## 1 column width 85 mm = 3.34 in
## 1.5 column width 114 mm = 4.49 in
## 2 columns width 174 mm = 6.85 in
## max height: 11 in

#########################################################################################################################
if (ref_sp == "human"){pdf_name="SupplementaryFigure20.pdf"; height=7}else{pdf_name="SupplementaryFigure21.pdf"; height=4.5}

pdf(paste(pathFigures, pdf_name, sep=""), height=height)

if (ref_sp == "human"){par(mfrow=c(4,1))}else{par(mfrow=c(2,1))}

#########################################################################################################################

## global synteny conservation, ENCODE

for (enh in enhancer.datasets[[ref_sp]]){
  m <- matrix(c(prop.obs.alldist[[enh]], prop.sim.alldist[[enh]]), nrow=2, byrow=T)
  
  par(mar=c(3.1, 4.1, 2.1, 1.1))
  
  if (enh == "FANTOM5"){ylim=c(50,105)}else{ylim = c(65,105)}
  if (ref_sp == "mouse"){ylim=c(65,100)}
  
  bar <-barplot(m, beside=T, space=c(0.25, 1.2), col=dataset.colors, border=dataset.colors, axes=F,  ylim = ylim, xpd=F)
  colnames(bar)<-species
  
  axis(side=2, las=2,  mgp=c(3, 0.75, 0), cex.axis=0.9)
  if (ref_sp == "mouse"){cex=0.85}else{cex=0.7}
  mtext("% pairs in \n conserved synteny", side=2, line=2, cex=cex)
  
  xax=apply(bar, 2, mean)
  axis(side=1, at=xax, mgp=c(3, 0.65, 0), labels=rep("", length(species)))
  mtext(species, side=1, at=xax, line=0.5, cex=0.75)
  
  mtext(paste(ref_sp, "vs."), side=1, line=0.5, at=min(xax)-diff(xax)[1]*1, cex=0.75)
  
  segments(bar[1,], ci.low.obs.alldist[[enh]], bar[1,], ci.high.obs.alldist[[enh]], col="black", lwd=1.1)
  
  segments(bar[2,], ci.low.sim.alldist[[enh]], bar[2,], ci.high.sim.alldist[[enh]], col="black", lwd=1.1)
  
  
  smallx=(bar[2,2]-bar[1,2])/10
  
  for (sp in species){
    this.pval=pvalues.alldist[[enh]][sp]
    
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
    
    ypos=max(prop.obs.alldist[[enh]][sp], prop.sim.alldist[[enh]][sp])+2
    
    segments(bar[1,sp]+smallx, ypos, bar[2,sp]-smallx, ypos)
    text(text, x=mean(as.numeric(bar[,sp])), y=ypos+2.5, xpd=NA, cex=1)
  }
  
  ## legend
  if (ref_sp == "mouse"){cexleg=0.85}else{cexleg=1}
  
  if (enh == "ENCODE"){
    legend("topright", legend = c("PCHi-C data", "simulated data"), fill=dataset.colors, border=dataset.colors, inset=c(0,-0.08), cex=cexleg, bty='n')}
  
  mtext(enh, side=3, cex=0.8, line=0.5)
}

#####################################################################################

dev.off()

#######################################################################################
