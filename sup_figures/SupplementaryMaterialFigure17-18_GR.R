#########################################################################################################################

objects=ls()

if(!"pathFigures"%in%objects){

 source("../main_figures/parameters.R")

 species.list=list()
 species.list[["human"]] <-c("macaque", "mouse", "rat",  "rabbit", "dog", "cow", "elephant", "opossum", "chicken")
 species.list[["mouse"]] <-c("rat", "rabbit", "macaque", "human", "dog", "cow", "elephant", "opossum", "chicken")
  
 enh="ENCODE" ## only one enhancer set
 
 load=T
}

#########################################################################################################################

if(load==TRUE){
  conserv_synteny_list=list()

  for (ref_sp in c("human", "mouse")){
    print(paste0("loading.. ", ref_sp))
    load(paste(pathFigures, "RData/data.synteny.conservation.",ref_sp,".RData", sep=""))
    conserv_synteny_list[[ref_sp]]=conserv_synteny
    rm("conserv_synteny")
  }
  
  load=FALSE
}

#########################################################################################################################

## 1 column width 85 mm = 3.34 in
## 1.5 column width 114 mm = 4.49 in
## 2 columns width 174 mm = 6.85 in
## max height: 11 in

## actual figure

for (ref_sp in c("human", "mouse")){
  species <-species.list[[ref_sp]]
  
  if (ref_sp == "human"){
    pdf_name="GenomeResearch_Figures/SupplementaryMaterialFigure17.pdf"
  } else{
    pdf_name="GenomeResearch_Figures/SupplementaryMaterialFigure18.pdf"
  }
  pdf(paste(pathFigures, pdf_name, sep=""), width=6.85, height=7.5)
  
  par(mfrow=c(3,3))
  par(mar=c(3.1, 4.1, 2.5, 1.1))
  
  ###################################################################################
  
  ## synteny conservation as a function of the distance between promoter and enhancer
  
  labels <- toupper(letters[1:length(species)])
  names(labels) <- species
  
  nbkept=25
  
  for (sp in species){
    conserv.obs <- conserv_synteny_list[[ref_sp]][[enh]][[sp]][["synt_obs"]]
    conserv.sim <- conserv_synteny_list[[ref_sp]][[enh]][[sp]][["synt_simul"]]
    
    prop.obs=tapply(conserv.obs$target_dist, conserv.obs$class_dist, function(x) 100*length(which(x<maxDistanceSyntenyTarget))/length(x))[1:nbkept]
    prop.sim=tapply(conserv.sim$target_dist, conserv.sim$class_dist, function(x) 100*length(which(x<maxDistanceSyntenyTarget))/length(x))[1:nbkept]
    
    ci.obs.low=tapply(conserv.obs$target_dist, conserv.obs$class_dist, function(x) 100*prop.test(length(which(x<maxDistanceSyntenyTarget)), length(x))$conf.int[1])[1:nbkept]
    ci.obs.high=tapply(conserv.obs$target_dist, conserv.obs$class_dist, function(x) 100*prop.test(length(which(x<maxDistanceSyntenyTarget)), length(x))$conf.int[2])[1:nbkept]
    
    ci.sim.low=tapply(conserv.sim$target_dist, conserv.sim$class_dist, function(x) 100*prop.test(length(which(x<maxDistanceSyntenyTarget)), length(x))$conf.int[1])[1:nbkept]
    ci.sim.high=tapply(conserv.sim$target_dist, conserv.sim$class_dist, function(x) 100*prop.test(length(which(x<maxDistanceSyntenyTarget)), length(x))$conf.int[2])[1:nbkept]
    
    par(mar=c(3.1, 4.1, 2.1, 1.1))
    
    ylim=range(c(ci.sim.low, ci.obs.low, ci.obs.high, ci.sim.high))
    
    smally=diff(ylim)/10
    ylim=ylim+c(-smally, smally)
    
    
    plot(1, type="n", xlim=c(0.5, nbkept+0.5), ylim=ylim, xlab="", ylab="", axes=F)
    
    xpos=1:nbkept
    
    ## lines(xpos, prop.obs, col=dataset.colors["Original"])
    points(xpos, prop.obs, col=dataset.colors["Original"], pch=20)
    segments(xpos, ci.obs.low, xpos, ci.obs.high, col=dataset.colors["Original"])
    
    ##lines(xpos, prop.sim, col=dataset.colors["Simulated"])
    points(xpos, prop.sim, col=dataset.colors["Simulated"], pch=20)
    segments(xpos, ci.sim.low, xpos, ci.sim.high, col=dataset.colors["Simulated"])
    
    axis(side=2, las=2,  mgp=c(3, 0.75, 0), cex.axis=0.9)
    
    xax=c(1, 5, 10, 15, 20, 25)
    axlab=as.character(c(0.05, 0.25, 0.5, 0.75, 1, 1.5))
    
    axis(side=1, mgp=c(3, 0.65, 0), at=xax, labels=axlab, cex.axis=0.95)
    
    mtext("% pairs in conserved synteny", side=2, line=2.5, cex=0.75)
    mtext("distance to promoters (Mb)", side=1, line=2, cex=0.75)

    
    ## legend
    if (labels[sp] == "A"){
      legend("bottomleft", legend = c("PCHi-C data", "simulated data"), col=dataset.colors, lty=1, cex=1.1, bty='n')
    }
        
    mtext(paste(ref_sp, "vs.", sp, sep=" "), side=3, cex=0.75, line=-0.5)

    mtext(labels[sp], side=3, cex=0.95, font=2, line=1, at=-7)
  }
  
  #####################################################################################
  
  dev.off()
  
  #######################################################################################
  
  
}

