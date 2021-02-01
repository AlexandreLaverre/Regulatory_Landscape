##########################################################################

library(Hmisc)

options(stringsAsFactors = FALSE)

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
  load(paste(pathFigures, "RData/data.sample.info.RData", sep=""))
  load(paste(pathFigures, "RData/data.enhancer.statistics.RData", sep=""))
  load(paste(pathFigures, "RData/data.contact.conservation.enhancers.RData", sep=""))

  # prepare data for each species
  
  chicago.dist <- list()
  chicago.dist.conf.low <- list()
  chicago.dist.conf.high <- list()
  
  cons.seq <- list()
  cons.seq.conf.low <- list()
  cons.seq.conf.high <- list()
  
  cons.contact <- list()
  cons.contact.conf.low <- list()
  cons.contact.conf.high <- list()
  
  cons.synteny <- list()
  cons.synteny.conf.low <- list()
  cons.synteny.conf.high <- list()

  nb_chicago_class = 10
  
  for(sp in c("human", "mouse")){
    ## prepare data

    tg=setdiff(c("human", "mouse"), sp)
    sampleinfo.tg=sampleinfo[[tg]]
      
    chicago.dist[[sp]] <- list()
    chicago.dist.conf.low[[sp]] <- list()
    chicago.dist.conf.high[[sp]] <- list()
    
    cons.seq[[sp]] <- list()
    cons.seq.conf.low[[sp]] <- list()
    cons.seq.conf.high[[sp]] <- list()
    
    cons.contact[[sp]] <- list()
    cons.contact.conf.low[[sp]] <- list()
    cons.contact.conf.high[[sp]] <- list()
    
    cons.synteny[[sp]] <- list()
    cons.synteny.conf.low[[sp]] <- list()
    cons.synteny.conf.high[[sp]] <- list()
    
    load(paste(pathFigures, "RData/data.sequence.conservation.stats.pcungapped.", sp, ".RData", sep=""))
    load(paste(pathFigures, "RData/data.synteny.conservation.", sp,".RData", sep=""))

    
    for (enh in enhancer.datasets[[sp]]){
      ## Calculate sequence conservation 
      enh_obs <- enhancer.statistics[[sp]][[enh]][["original"]]
      enh_obs$class_score <- cut2(enh_obs$median_score, g=nb_chicago_class)
      
      enh_align <- list_align_enh[[enh]][["enh_align_obs"]]
      rownames(enh_align) <- enh_align$ID
      common=intersect(rownames(enh_obs), rownames(enh_align))
      
      enh_align = enh_align[common,]
      enh_align$class_score = enh_obs[common,]$class_score
      
      cons.seq[[sp]][[enh]] <- tapply(100*enh_align[[tg]], as.factor(enh_align$class_score), function(x) mean(x, na.rm=T))
      cons.seq.conf.low[[sp]][[enh]] <- tapply(100*enh_align[[tg]], as.factor(enh_align$class_score), function(x) t.test(x)[["conf.int"]][1])
      cons.seq.conf.high[[sp]][[enh]] <- tapply(100*enh_align[[tg]], as.factor(enh_align$class_score), function(x) t.test(x)[["conf.int"]][2])
      
      ## Calculate relation with distance to promoter
      contact_obs <- contact.conservation[[paste0(sp, "2", tg)]][[enh]][["obsobs"]] 
      
      contact_obs$class_dist <- cut(contact_obs$origin_dist, breaks=seq(from=minDistance, to=maxDistance, by=50000), include.lowest = T)
      
      chicago.dist[[sp]][[enh]] <- tapply(contact_obs$origin_med_score, as.factor(contact_obs$class_dist), function(x) mean(x, na.rm=T))
      chicago.dist.conf.low[[sp]][[enh]] <- tapply(contact_obs$origin_med_score, as.factor(contact_obs$class_dist), function(x) t.test(x)[["conf.int"]][1])
      chicago.dist.conf.high[[sp]][[enh]] <- tapply(contact_obs$origin_med_score, as.factor(contact_obs$class_dist), function(x) t.test(x)[["conf.int"]][2])
      
      ## Calculate contact conservation proportion
      contact_obs$cons=apply(contact_obs[,sampleinfo.tg$Sample.ID], 1, function(x) any(x>0))
      
      contact_obs$class_score <- cut2(contact_obs$origin_med_score, g=nb_chicago_class)
      
      cons.contact[[sp]][[enh]] <- sapply(levels(contact_obs$class_score), function(x)
                                    (nrow(contact_obs[which(contact_obs$cons == T & contact_obs$class_score == x),])/nrow(contact_obs[which(contact_obs$class_score == x ),]))*100)
      
      cons.contact.conf.low[[sp]][[enh]] <- sapply(levels(contact_obs$class_score), function(x)
                                             prop.test(nrow(contact_obs[which(contact_obs$cons == T & contact_obs$class_score == x),]), nrow(contact_obs[which(contact_obs$class_score == x ),]))$conf.int[1]*100)
      
      cons.contact.conf.high[[sp]][[enh]] <- sapply(levels(contact_obs$class_score), function(x)
                                            prop.test(nrow(contact_obs[which(contact_obs$cons == T & contact_obs$class_score == x),]), nrow(contact_obs[which(contact_obs$class_score == x ),]))$conf.int[2]*100)
      
    ## Calculate synteny conservation 
      synt_obs <- conserv_synteny[[enh]][[tg]][["synt_obs"]]
      rownames(synt_obs) <- paste0(synt_obs$origin_gene, "-", synt_obs$origin_enh)
      rownames(contact_obs) <- paste0(contact_obs$origin_gene, "-", contact_obs$origin_enh)
      
      common=intersect(rownames(synt_obs), rownames(contact_obs))
      synt_obs <- synt_obs[common,]
      synt_obs$class_score = contact_obs[common,]$class_score
      
      cons.synteny[[sp]][[enh]] <- sapply(levels(synt_obs$class_score), function(x)
                                    (nrow(synt_obs[which(synt_obs$target_dist <= maxDistanceSyntenyTarget & synt_obs$class_score == x),])/nrow(synt_obs[which(synt_obs$class_score == x ),]))*100)
      
      cons.synteny.conf.low[[sp]][[enh]] <- sapply(levels(contact_obs$class_score), function(x)
                                             prop.test(nrow(synt_obs[which(synt_obs$target_dist <= maxDistanceSyntenyTarget & synt_obs$class_score == x),]), nrow(synt_obs[which(synt_obs$class_score == x ),]))$conf.int[1]*100)
      
      cons.synteny.conf.high[[sp]][[enh]] <- sapply(levels(contact_obs$class_score), function(x)
                                              prop.test(nrow(synt_obs[which(synt_obs$target_dist <= maxDistanceSyntenyTarget & synt_obs$class_score == x),]), nrow(synt_obs[which(synt_obs$class_score == x ),]))$conf.int[2]*100)
      
    }
    
  }
  
  load=FALSE
}

#####################################################################################
## make the actual figure

for(sp in c("human", "mouse")){
   
  if (sp == "human"){
    pdf.name = "SupplementaryFigure25.pdf"
  }else{
    pdf.name = "SupplementaryFigure26.pdf"
  }
  
  print(pdf.name)
  
  tg=setdiff(c("human", "mouse"), sp)
  sampleinfo.tg=sampleinfo[[tg]]
  
  ###################################################################################
  
  pdf(paste(pathFigures, pdf.name, sep=""), width=6.85, height=5.5)
  
  par(mai = c(0.5, 0.5, 0.5, 0.2)) # bottom, left, top, right
  mtext.CEX = 0.75
  
  m=matrix(rep(NA, 2*10), nrow=2)
  m[1,]=c(rep(1,5), rep(2,5))
  m[2,]=c(rep(3,5), rep(4,5))
  layout(m)
  
############################## CHICAGO score according to dist ###############################
  
  if (sp == "human"){
    ylim=c(6, 12)
  }else{
    ylim=c(6, 9)
  }
  
  plot(1, type="n", xlim=c(0.5, 40.5), ylim=ylim, xlab="", ylab="", axes=F)
  xpos=1:39
  
  for(enh in enhancer.datasets[[sp]]){
    lines(xpos, chicago.dist[[sp]][[enh]][xpos], col=col.enhancers[enh])
    segments(xpos, chicago.dist.conf.low[[sp]][[enh]][xpos], xpos, chicago.dist.conf.high[[sp]][[enh]][xpos], col=col.enhancers[enh])
  }
  
  xax=c(1, 10,  20, 30, 40)
  axlab=as.character(c(0.05,  0.5,  1, 1.5, 2))
  axis(side=1, mgp=c(3, 0.65, 0), at=xax, labels=axlab)
  mtext("distance to promoters (Mb)", side=1, line=2, cex=mtext.CEX)
  
  axis(side=2, mgp=c(3, 0.75, 0))
  mtext("CHICAGO score", side=2, line=2.5, cex=mtext.CEX)
  
  legend("topleft", legend=label.enhancers[enhancer.datasets[[sp]]], lty=1, 
         col=col.enhancers[enhancer.datasets[[sp]]], bty="n", cex=1, seg.len=1)

  ## plot label
  mtext("a", side=3, font=2, at=-7.2, line=2.5)
  
  
  ############################## Sequence conservation ###############################
  xpos=seq(1, nb_chicago_class, 1)
  names(xpos) = 1:nb_chicago_class
  smallx=c(-0.15, -0.075, 0.075, 0.15)
  names(smallx)=enhancer.datasets[[sp]]
  
  if (sp == "human"){
    ylim=c(30, 60)
  }else{
    ylim=c(40, 60)
  }
  
  xlim=c(0.5, nb_chicago_class+0.5)
  
  plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=ylim, xaxs="i", yaxs="i")
  
  for(enh in enhancer.datasets[[sp]]){
    for(row in 1:length(cons.seq[[sp]][[enh]])){
      
      x=xpos[row]+smallx[enh]
      points(x, cons.seq[[sp]][[enh]][row], pch=20, col=col.enhancers[enh], cex=1.1)
      segments(x, cons.seq.conf.low[[sp]][[enh]][row], x, cons.seq.conf.high[[sp]][[enh]][row], col=col.enhancers[enh])
    }
  }
  
  abline(v=xpos[1:nb_chicago_class-1]+0.5, lty=3, col="gray40")
  axis(side=1, at=xpos, mgp=c(3, 0.5, 0), labels=rep("", nb_chicago_class))
  mtext(1:nb_chicago_class, at=xpos, side=1, line=1, cex=mtext.CEX)
  mtext("CHICAGO score decile", side=1, line=2.5, cex=mtext.CEX)
  
  axis(side=2, mgp=c(3, 0.75, 0))
  mtext(paste0("% aligned sequence in ", tg), side=2, line=2.5, cex=mtext.CEX)

  ## plot label
  mtext("b", side=3, font=2, at=-0.9, line=2.5)
  
  ############################## Synteny conservation ###############################
  if (sp == "human"){
    ylim=c(93, 100)
  }else{
    ylim=c(92, 100)
  }
  
  plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=ylim, xaxs="i", yaxs="i")
  
  for(enh in enhancer.datasets[[sp]]){
    for(row in 1:length(cons.synteny[[sp]][[enh]])){
      
      x=xpos[row]+smallx[enh]
      points(x, cons.synteny[[sp]][[enh]][row], pch=20, col=col.enhancers[enh], cex=1.1)
      segments(x, cons.synteny.conf.low[[sp]][[enh]][row], x, cons.synteny.conf.high[[sp]][[enh]][row], col=col.enhancers[enh])
    }
  }
  
  abline(v=xpos[1:nb_chicago_class-1]+0.5, lty=3, col="gray40")
  axis(side=1, at=xpos, mgp=c(3, 0.5, 0), labels=rep("", nb_chicago_class))
  mtext(1:nb_chicago_class, at=xpos, side=1, line=1, cex=mtext.CEX)
  mtext("CHICAGO score decile", side=1, line=2.5, cex=mtext.CEX)
  
  axis(side=2, mgp=c(3, 0.75, 0))
  mtext(paste0("% pairs in conserved synteny in ", tg), side=2, line=2.5, cex=mtext.CEX)
  
  ## plot label
  mtext("c", side=3, font=2, at=-0.9, line=2.5)
  
############################## Contact conservation ###############################
  
  if (sp == "human"){
    ylim=c(10, 60)
  }else{
    ylim=c(20, 55)
  }
    
  plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=ylim, xaxs="i", yaxs="i")
  
  for(enh in enhancer.datasets[[sp]]){
    for(row in 1:length(cons.contact[[sp]][[enh]])){
      
      x=xpos[row]+smallx[enh]
      points(x, cons.contact[[sp]][[enh]][row], pch=20, col=col.enhancers[enh], cex=1.1)
      segments(x, cons.contact.conf.low[[sp]][[enh]][row], x, cons.contact.conf.high[[sp]][[enh]][row], col=col.enhancers[enh])
    }
  }
  
  abline(v=xpos[1:nb_chicago_class-1]+0.5, lty=3, col="gray40")
  axis(side=1, at=xpos, mgp=c(3, 0.5, 0), labels=rep("", nb_chicago_class))
  mtext(1:nb_chicago_class, at=xpos, side=1, line=1, cex=mtext.CEX)
  mtext("CHICAGO score decile", side=1, line=2.5, cex=mtext.CEX)
  
  axis(side=2, mgp=c(3, 0.75, 0))
  mtext(paste0("% pairs in conserved contact in ", tg), side=2, line=2.5, cex=mtext.CEX)

  
  ## plot label
  mtext("d", side=3, font=2, at=-0.9, line=2.5)
  
  #####################################################################################
  
  dev.off()
  
  #####################################################################################
  
}

