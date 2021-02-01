######################################################################################################################

objects=ls()

if(!"pathScripts"%in%objects){
  load=T
  prepare=T
  source("parameters.R")
}

##############################################################################

if(load){
  sp="human"
  
  load(paste(pathFigures, "RData/data.gene.annotations.RData", sep=""))
  load(paste(pathFigures, "RData/data.", sp, ".CM2019.SomaticOrgans.expdiv.RData", sep=""))
  load(paste(pathFigures, "RData/data.", sp, ".regland.conservation.RData", sep=""))
  
  if(sp == "human"){
    sp_name="Human"
  } else{
    sp_name="Mouse"
  }
  
  load=FALSE
  
}

#############################################################################################################
######################## Part 1 : Complexity and gene expression structure ##################################

xpos=seq(1, 5, 1)
smallx=c(-0.15, -0.075, 0.075, 0.15)
names(smallx)=enhancer.datasets[[sp]]
xlim=c(0.5, 5.5)
cex.mtext = 0.75
syn.datasets=c("Original", "Simulated")
names(syn.datasets)=c("obs", "sim")

smallxdata=c(-0.075, 0.075)
names(smallxdata)=c("obs", "sim")

smallxdist=c(-0.2, 0, 0.2)
names(smallxdist)=c("all", "25kb - 100kb", "100kb - 500kb")

distances=c("all", "25kb - 100kb", "100kb - 500kb")
col.distances=c("black", "steelblue", "indianred")
names(col.distances)=distances

CMPlot <- function(var, plot.nb){
  if (var == "CorrelationSpearman"){
    ylim=c(0.55, 0.65)
  } else if (var == paste0("Tau", sp_name)){
    ylim=c(0.6, 0.77)
  } else if (var == "MeanRPKM"){
    ylim=c(8, 11.25)
  } else if (var == "CorrectedSpearman"){
    ylim=c(0.02, 0.1)
  }
  
  plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=ylim, xaxs="i", yaxs="i")

  type="obs"
  
  for(dist in distances){
    for(enh in c("ENCODE")){
      regland = genes.conservation[[enh]][[type]][[dist]]
      genes = intersect(rownames(regland), rownames(expdiv))
      regland = regland[genes,]

  
      for(class in levels(regland$class_nb_contact)){
        this.genes=rownames(regland[which(regland$class_nb_contact == class),])
        
        xpos=seq(1, length(levels(regland$class_nb_contact)), 1)
        names(xpos) = levels(regland$class_nb_contact)
 
        x=xpos[class]+smallxdist[dist]
        
        b=boxplot(expdiv[this.genes, var], plot=FALSE)
        med=median(expdiv[this.genes, var])
        ci=as.numeric(b$conf)

        points(x, med, pch=20, col=col.distances[dist], cex=1.1, type="b")
        segments(x, ci[1], x, ci[2],  col=col.distances[dist])
      }
    }
  }
  
  abline(v=xpos[1:4]+0.5, lty=3, col="gray40")
  
  axis(side=2, mgp=c(3, 0.75, 0), cex.axis=1, las=2)
  mtext(names_MeasuresCM[plot.nb], side=2, line=2.95, cex=cex.mtext)
  
  ## plot label
  mtext(letters[plot.nb], side=3, line=1.5, at=-0.85, font=2, cex=1.2)
  
  axis(side=1, cex.axis=1, mgp=c(3, 0.75, 0), at=xpos, labels=rep("",5))
  mtext(c("low", "", "medium", "", "high"), at=xpos, side=1, line=0.75, cex=0.75) 
  mtext("number of contacts", side=1, line=2.5, cex=cex.mtext)
  
 
}

################################################################################################################################
############################## PART2 :  Cardoso-Moreira plot function ##########################################################

plot_profiles <- function(class_conserv, distances, xlab, xnames){
  type="obs"
  enh="ENCODE"
  
  smallx=c(-0.15, -0.075, 0.075, 0.15)
  names(smallx)=enhancer.datasets[[sp]]
  
  if (Measure == "corrected"){
    DivergenceMeasure = "CorrectedSpearman"
    ylab="expression conservation"
  
  }else{
    DivergenceMeasure = "CorrelationSpearman"
    ylab="expression conservation"
  }

  ## define ylim first

  all.values=c()

  for (dist in distances){
    regland = genes.conservation[[enh]][[type]][[dist]]
    genes = intersect(rownames(regland), rownames(expdiv))
    regland = regland[genes,]
    
    for (class in levels(regland[[class_conserv]])){
      this.genes=rownames(regland[which(regland[[class_conserv]] == class),])
      
      b=boxplot(expdiv[this.genes, DivergenceMeasure], plot=F)
      med=median(expdiv[this.genes, DivergenceMeasure], na.rm=T)
      ci=as.numeric(b$conf)
      all.values=c(all.values, ci)
    }
  }

  ylim=range(all.values)
  smally=diff(ylim)/5
  ylim=ylim+c(-smally,smally)
  
  ## now do the actual plot
 
  xmax=length(xnames)
  
  xlim=c(0.5, xmax+0.5)
  
  plot(1, type="n", xlab="", ylab="", axes=F, xlim=xlim, ylim=ylim, xaxs="i", yaxs="i")
  
  for (dist in distances){
    regland = genes.conservation[[enh]][[type]][[dist]]
    genes = intersect(rownames(regland), rownames(expdiv))
    regland = regland[genes,]
    
    for (class in levels(regland[[class_conserv]])){
      this.genes=rownames(regland[which(regland[[class_conserv]] == class),])
      
      b=boxplot(expdiv[this.genes, DivergenceMeasure], plot=F)
      med=median(expdiv[this.genes, DivergenceMeasure], na.rm=T)
      ci=as.numeric(b$conf)
      
      xpos=seq(1,  length(levels(regland[[class_conserv]])), 1)
      names(xpos) = levels(regland[[class_conserv]])
      
      x=xpos[class]+smallxdist[dist]
      
      points(x, med, pch=20, col=col.distances[dist], cex=1.1)
      segments(x, ci[1], x, ci[2], col=col.distances[dist])
      
    }
  }
  
  
  abline(v=xpos[1:xmax-1]+0.5, lty=3, col="gray40")
  axis(side=1, at=xpos, mgp=c(3, 0.75, 0), labels=rep("", length(levels(regland[[class_conserv]]))), cex.axis=cex.mtext)
  
  if(class_conserv=="class_align_score"){
    mtext(c("low", "", "medium", "", "high"), at=xpos, side=1, line=0.75, cex=cex.mtext)
  } else{
   mtext(xnames, at=xpos, side=1, line=0.75, cex=cex.mtext)
 }
   mtext(xlab, side=1, line=2.5, cex=cex.mtext)
  
  axis(side=2, mgp=c(3, 0.75, 0), cex.axis=1, las=2)
  mtext(ylab, side=2, line=2.95, cex=cex.mtext)
    
}

################################################################################################################################
################################################### Plot Figure 6 ##############################################################

Measure = "corrected"

pdf(file=paste(pathFigures, "Figure6.pdf", sep=""), width = 6.85,  height=6)

m=matrix(rep(NA,2*9), nrow=2)
m[1,]=c(rep(c(1,2,3), each=3))
m[2,]=c(rep(4,3), rep(5,2), rep(6,2), rep(7,2))

layout(m)

par(mar = c(5.1, 4.5, 2.75, 0.5))

######################## Part 1 : Complexity ########################

MeasuresCM <- c("MeanRPKM", paste0("Tau", sp_name),  "CorrectedSpearman")
names_MeasuresCM <- c("mean expression level (RPKM)", "expression specificity", "expression conservation")

for (measure in 1:length(MeasuresCM)){
  CMPlot(MeasuresCM[measure], measure)
}

######################## Part 2 : Gene expression profile evolution ################

plot_profiles("class_align_score", distances,  "enhancer sequence conservation", 1:5)
mtext("d", side=3, at=-0.85, font=2, cex=1.2, line=1)

plot_profiles("class_cons_synt", distances,  "% conserved synteny", c("<100%", "100%"))
mtext("e", side=3, at=-0.5, font=2, cex=1.2, line=1)

plot_profiles("class_cons_cont", distances,  "% conserved contacts", c("<50%",">50%"))
mtext("f", side=3, at=-0.5, font=2, cex=1.2, line=1)

################################################################################################################################
## empty plot for the legend
par(mar=c(0, 0.1, 0, 0.3)) 
plot.new()

legend("topleft", col=c(col.distances[1], "white", col.distances[2], "white", col.distances[3]), legend = c("all", "", "short range\n(25 kb - 100 kb)","", "medium range\n(100 kb - 500 kb)"), box.col="white", bg="white", pch=20, cex=1.1, inset=c(0.01,0.1), xpd=T)

## legend("topleft", col=c(col.distances[1], "white", col.distances[2]), legend = c("short range\n(25 kb - 100 kb)","", "medium range\n(100 kb - 500 kb)"), box.col="white", bg="white", pch=20, cex=1.1, inset=c(0.01,0.1), xpd=T)

################################################################################################################################

dev.off()

################################################################################################################################
