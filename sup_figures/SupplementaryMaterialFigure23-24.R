#########################################################################################################################

objects=ls()

if(!"pathFigures"%in%objects){

 source("../main_figures/parameters.R")

 species.list=list()
 species.list[["human"]] <-c("macaque", "mouse", "rat",  "rabbit", "dog", "cow", "elephant", "opossum", "chicken")
 species.list[["mouse"]] <-c("rat", "rabbit", "macaque", "human", "dog", "cow", "elephant", "opossum", "chicken")
 
 load=T
}

#########################################################################################################################

if(load){
  prop.obs.alldist <- list()
  prop.sim.alldist <- list()
  ci.low.obs.alldist <- list()
  ci.low.sim.alldist <- list()
  ci.high.obs.alldist <- list()
  ci.high.sim.alldist <- list()
  pvalues.alldist <- list()
  
  for (ref_sp in c("human", "mouse")){
    print(paste0('loading... ', ref_sp))
    
    load(paste(pathFigures, "RData/data.synteny.conservation.",ref_sp,".RData", sep=""))
  
    species=species.list[[ref_sp]]
    
    prop.obs.alldist[[ref_sp]] <- list()
    prop.sim.alldist[[ref_sp]] <- list()
    ci.low.obs.alldist[[ref_sp]] <- list()
    ci.low.sim.alldist[[ref_sp]] <- list()
    ci.high.obs.alldist[[ref_sp]] <- list()
    ci.high.sim.alldist[[ref_sp]] <- list()
    pvalues.alldist[[ref_sp]] <- list()
       
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
        
        pvalues.alldist[[ref_sp]][[enh]] <- c(pvalues.alldist[[ref_sp]][[enh]], pval)
        
        prop.obs.alldist[[ref_sp]][[enh]] <- c(prop.obs.alldist[[ref_sp]][[enh]], 100*nb_cons_synt_obs/nb_tot_synt_obs)
        prop.sim.alldist[[ref_sp]][[enh]] <- c(prop.sim.alldist[[ref_sp]][[enh]], 100*nb_cons_synt_simul/nb_tot_synt_simul)
        
        prop.test.obs <- prop.test(x = nb_cons_synt_obs, n=nb_tot_synt_obs, p=0.5)
        prop.test.simul <- prop.test(x = nb_cons_synt_simul, n=nb_tot_synt_simul, p=0.5)
        
        ci.low.obs.alldist[[ref_sp]][[enh]] <- c(ci.low.obs.alldist[[ref_sp]][[enh]], 100*prop.test.obs$conf.int[1])
        ci.high.obs.alldist[[ref_sp]][[enh]] <- c(ci.high.obs.alldist[[ref_sp]][[enh]], 100*prop.test.obs$conf.int[2])
        
        ci.low.sim.alldist[[ref_sp]][[enh]] <- c(ci.low.sim.alldist[[ref_sp]][[enh]], 100*prop.test.simul$conf.int[1])
        ci.high.sim.alldist[[ref_sp]][[enh]] <- c(ci.high.sim.alldist[[ref_sp]][[enh]], 100*prop.test.simul$conf.int[2])
      }
      
      names(prop.obs.alldist[[ref_sp]][[enh]]) <- species
      names(prop.sim.alldist[[ref_sp]][[enh]]) <- species
      
      names(ci.low.obs.alldist[[ref_sp]][[enh]]) <- species
      names(ci.high.obs.alldist[[ref_sp]][[enh]]) <- species
      
      names(ci.low.sim.alldist[[ref_sp]][[enh]]) <- species
      names(ci.high.sim.alldist[[ref_sp]][[enh]]) <- species
      
      names(pvalues.alldist[[ref_sp]][[enh]]) <- species
    }
  }
  
  ## we do this just once
  load=FALSE
}

#########################################################################################################################

## 1 column width 85 mm = 3.34 in
## 1.5 column width 114 mm = 4.49 in
## 2 columns width 174 mm = 6.85 in
## max height: 11 in

for (ref_sp in c("human", "mouse")){
 
  if (ref_sp == "human"){
    pdf_name="SupplementaryMaterialFigure23.pdf"
    height=7
    cex.axis=0.9
    cex.lab=0.85
    cexleg=1
    cex.text=1
    xpart=7.5
  } else{
    pdf_name="SupplementaryMaterialFigure24.pdf"
    height=4.5
    cex.axis=0.7
    cex.lab=0.7
    cexleg=0.7
    cex.text=0.7
    xpart=6.5
  }
  
  pdf(paste(pathFigures, pdf_name, sep=""), height=height)
  
  if (ref_sp == "human"){
    par(mfrow=c(4,1))
     par(mar=c(3.1, 4.5, 2.1, 1.1))
  }else{
    par(mfrow=c(2,1))
     par(mar=c(3.1, 3.7, 2.1, 1.1))
  }
  
#########################################################################################################################

  species <- species.list[[ref_sp]]
  
  labels=letters[1:length(enhancer.datasets[[ref_sp]])]
  names(labels)=enhancer.datasets[[ref_sp]]

  for (enh in enhancer.datasets[[ref_sp]]){
    m <- matrix(c(prop.obs.alldist[[ref_sp]][[enh]], prop.sim.alldist[[ref_sp]][[enh]]), nrow=2, byrow=T)
    
    ylim=c(50,105)
      
    bar <-barplot(m, beside=T, space=c(0.25, 1.2), col=dataset.colors, border=dataset.colors, axes=F,  ylim = ylim, xpd=F)
    colnames(bar)<-species
    
    axis(side=2, las=2,  mgp=c(3, 0.75, 0), cex.axis=cex.axis)
   
    mtext("% pairs in \n conserved synteny", side=2, line=2, cex=cex.lab)
    
    xax=apply(bar, 2, mean)
    axis(side=1, at=xax, mgp=c(3, 0.65, 0), labels=rep("", length(species)))
    mtext(species, side=1, at=xax, line=0.5, cex=0.75)
    
    mtext(paste(ref_sp, "vs."), side=1, line=0.5, at=min(xax)-diff(xax)[1]*1, cex=0.75)
    
    segments(bar[1,], ci.low.obs.alldist[[ref_sp]][[enh]], bar[1,], ci.high.obs.alldist[[ref_sp]][[enh]], col="black", lwd=1.1)
    segments(bar[2,], ci.low.sim.alldist[[ref_sp]][[enh]], bar[2,], ci.high.sim.alldist[[ref_sp]][[enh]], col="black", lwd=1.1)
    
    smallx=(bar[2,2]-bar[1,2])/10
    
    for (sp in species){
      this.pval=pvalues.alldist[[ref_sp]][[enh]][sp]
      
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
      
      ypos=max(prop.obs.alldist[[ref_sp]][[enh]][sp], prop.sim.alldist[[ref_sp]][[enh]][sp])+2
      
      segments(bar[1,sp]+smallx, ypos, bar[2,sp]-smallx, ypos)
      text(text, x=mean(as.numeric(bar[,sp])), y=ypos+2.5, xpd=NA, cex=cex.text)
    }
    
    ## legend
 
    
    if (enh == "ENCODE"){
      legend("topright", legend = c("PCHi-C data", "simulated data"), fill=dataset.colors, border=dataset.colors, inset=c(0,-0.08), cex=cexleg, bty='n')
    }
    
    mtext(enh.syn[enh], side=3, cex=0.8, line=0.5)

    xlim=range(as.numeric(bar))
    
    mtext(labels[enh], side=3, font=2, line=0.95, at=xlim[1]-diff(xlim)/xpart)
  }
  
  dev.off()
  
}

#########################################################################################################################

