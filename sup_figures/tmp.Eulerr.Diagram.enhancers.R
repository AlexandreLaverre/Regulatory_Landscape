###########################################################################################
## if it's the first time we run this figure, we load and prepare data
library(eulerr)
options(stringsAsFactors = FALSE)

objects=ls()

if(!"pathScripts"%in%objects){
  prepare=T
  source("parameters.R") ## paths are defined based on the user name
  path <- "/home/laverre/Documents/Regulatory_Landscape/data/"
  species <- c("human", "mouse")
}

###########################################################################################

if (prepare){
  overlap <- list()
  
  for (ref_sp in species){
    overlap[[ref_sp]] <- read.table(paste(path, ref_sp, "/potential_enhancers/enhancers_coordinates/common_enhancers_overlap", sep=""), header=T, row.names=1)
    
    overlap[[ref_sp]][!is.na(overlap[[ref_sp]])] <- TRUE
    overlap[[ref_sp]][is.na(overlap[[ref_sp]])] <- FALSE
    
    colnames(overlap[[ref_sp]]) <- enhancer.datasets[[ref_sp]]
    
    overlap[[ref_sp]] <- apply(overlap[[ref_sp]], 2, function(x) as.logical(x))
    overlap[[ref_sp]]$ENCODE <- as.logical(overlap[[ref_sp]]$ENCODE)

  }
  prepare = FALSE
}

###########################################################################################
########## output ##########
pdf(paste(pathFigures, "/SupplementaryFigureX_overlap_between_enhancers.pdf", sep=""), width=7, height=7)

par(mfrow=c(1,2))
par(mai = c(0.5, 0.6, 0.3, 0.4)) # bottom, left, top, right

for (ref_sp in species){
  plot(euler(overlap[[ref_sp]]), quantities = TRUE, fill = col.enhancers)
}

dev.off()
