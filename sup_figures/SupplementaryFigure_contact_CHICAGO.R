################################################################################################################################################
library(data.table)
options(stringsAsFactors = FALSE)

#setwd("/home/laverre/Data/Regulatory_landscape/scripts/main_figures")
#########################################################################

## if it's the first time we run this figure, we load and prepare data

objects=ls()

if(!"pathScripts"%in%objects){
  load=T
  prepare=T
  source("../main_figures/parameters.R")
}

##################################################################

if(load){
  sp="human"
  
  # Choose genes within : all ; dvpt ; other
  selected_genes = "all" 
  if (selected_genes == "dvpt"){`%get%` <- `%in%`}else{`%get%` <- Negate(`%in%`)}
  
  minDistance=25e3
  maxDistance=2e6
  
  load(paste(pathFigures, "RData/data.gene.enhancer.contacts.RData", sep=""))
  load(paste(pathFigures, "RData/data.gene.enhancer.contacts.conservation.RData", sep=""))
  load(paste(pathFigures, "RData/data.enhancer.statistics.RData", sep=""))
  
  path_overlap=paste(pathFinalData, "SupplementaryDataset7/", sp, "/sequence_conservation/enhancers/", sep="")
  
  if(selected_genes != "all"){
    dev_gene <- read.table(paste(pathFinalData, "SupplementaryDataset3/genes/", sp, "_dvpt_process_genes.txt", sep=""), header=T, sep="\t")
  }
  
  load=FALSE
}

##################################################################

for (enh in enhancer.datasets[[sp]]){
  message(enh)
  all_obs <- gene.enhancer.contacts[[sp]][[enh]][["real"]] 
  contact_obs <- contacts.conservation[[sp]][[enh]][["real"]] 

  
  # Select only contacts with no duplicated enhancers
  obs_stats <- enhancer.statistics[[sp]][[enh]][["real"]]
  obs_stats$enh <-  do.call(paste,c(obs_stats[c("chr","start","end")],sep=":"))
  obs_stats <- obs_stats[which(obs_stats$BLAT_match < 2),]

  all_obs <- all_obs[which(all_obs$enhancer %in% obs_stats$enh),]
  contact_obs <- contact_obs[which(contact_obs$origin_enh %in% obs_stats$enh),]

  
  # Sample classes
  all_obs$class_score <- cut(all_obs$median_score, breaks=c(5, 6, 7, 8, 9,10, max(contact_obs$origin_med_score)), include.lowest = T)
  contact_obs$class_score <- cut(contact_obs$origin_med_score, breaks=c(5, 6, 7, 8, 9,10, max(contact_obs$origin_med_score)), include.lowest = T)

  # Calculate proportion
  conserv_obs <- data.frame(result = sapply(levels(all_obs$class_score), function(x)
    (nrow(contact_obs[which(contact_obs$class_score == x),])/(nrow(all_obs[which(all_obs$class_score == x ),])))*100))
  
  conserv_obs$conflow <- sapply(levels(all_obs$class_score), function(x)  
    (prop.test(x = nrow(contact_obs[which(contact_obs$class_score == x ),]),
               n=nrow(all_obs[which(all_obs$class_score == x ),]), p=0.5)$conf.int[1])*100)
  
  conserv_obs$confup <- sapply(levels(all_obs$class_score), function(x)  
    (prop.test(x = nrow(contact_obs[which(contact_obs$class_score == x ),]),
               n=nrow(all_obs[which(all_obs$class_score == x ),]), p=0.5)$conf.int[2])*100)
  
  
  conserv_obs$count <- sapply(levels(all_obs$class_score), function(x) paste("n=", nrow(all_obs[which(all_obs$class_score == x ),]), sep=""))
  
  assign(paste("conserv_sample", enh, sep = "_"), conserv_obs)
  
}

########################################################################################
pdf(paste(pathFigures, "SupplementaryFigureX_conservcontact_CHICAGOscore.pdf", sep=""), width=7, height=5)

ylim=c(0, 20)
xlim=c(0.5, 6.5)

xpos=seq(1, 6, 1)
names(xpos) = 1:6

smallx=c(-0.15, -0.075, 0.075, 0.15)
names(smallx)=enhancer.datasets[[sp]]

plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=ylim, xaxs="i", yaxs="i")

for(enh in enhancer.datasets[[sp]]){
  conserv = get(paste("conserv_sample", enh, sep="_"))
  
  for(row in 1:nrow(conserv)){

    x=xpos[row]+smallx[enh]
    points(x, conserv[row,'result'], pch=20, col=col.enhancers[enh], cex=1.1)
    segments(x, conserv[row,'conflow'], x, conserv[row,'confup'], col=col.enhancers[enh])
  }
}

abline(v=xpos[1:5]+0.5, lty=3, col="gray40")
axis(side=1, at=xpos, mgp=c(3, 0.5, 0), labels=rep("", 6))
mtext(levels(all_obs$class_score), at=xpos, side=1, line=1)
mtext("CHICAGO Score", side=1, line=2.5)

axis(side=2, mgp=c(3, 0.75, 0))
mtext("Conserved contact (%)", side=2, line=2.5)

legend("topleft", legend=label.enhancers[enhancer.datasets[[sp]]], lty=1, 
       col=col.enhancers[enhancer.datasets[[sp]]], bty="n", cex=1, seg.len=1)

dev.off()
