
###########################################################################################
objects=ls()

if(!"pathScripts"%in%objects){
  source("../main_figures/parameters.R") ## paths are defined based on the user name
  nb_draw = 100
  cumul_int = list()
}

###########################################################################################

for (sp in c("human", "mouse")){
  for (data in c("observed", "simulated")){
    print(paste("loading...", sp, "in", data, "data", sep=" "))
    
    if (data == "observed"){
      setwd(paste(pathFinalData, "SupplementaryDataset1", sp, "interactions_samples", sep="/"))
    } else{
      setwd(paste(pathFinalData, "SupplementaryDataset2", sp, "interactions_samples", sep="/"))
    }
    
    # Load all interactions
    all.files=system("ls", intern=T)
    all.int = list()
    for (file in all.files){
      int = read.table(file, h=T)
      int = int[which(int$baited_frag == "unbaited" & int$dist <= maxDistance & int$dist >= minDistance),]
      int$ID <-  do.call(paste,c(int[c("bait_chr","bait_start","bait_end","chr","start","end")],sep=":"))
      all.int[[file]] = int$ID
    }
    
    
    # Count cumulative interactions
    cumul_int[[sp]][[data]] = matrix(nrow=length(all.files), ncol=nb_draw)
    
    for(i in 1:nb_draw){
      sampled.files  = sample(all.files, length(all.files), replace=F)
      interactions = c()
      n_sample = 1
      
      for (sample in sampled.files){
        interactions = c(interactions, all.int[[sample]][which(!all.int[[sample]] %in% interactions)])
        
        cumul_int[[sp]][[data]][n_sample,i] = length(interactions)
        
        n_sample = n_sample+1
      }
    }
  }
}

###########################################################################################

### Output
save(cumul_int, file = paste(pathFigures, "RData/data.samples.cumulative.interactions.RData", sep=""))

###########################################################################################
