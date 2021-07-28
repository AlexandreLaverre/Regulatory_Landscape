###########################################################################

library(data.table)

source("parameters.R") ## paths are defined based on the user name

pathPhyloP="/home/laverre/Regulatory_landscape/result/phyloP/"

###########################################################################

load(paste(pathFigures, "RData/data.enhancer.statistics.RData", sep=""))
load(paste(pathFigures, "RData/data.fragment.statistics.RData", sep=""))

###########################################################################

minAlnLength=10

###########################################################################

for(ref in c("human", "mouse")){
  if (ref == "human"){way="100way"}else{way="60way"}
  
  frag.obs=fragment.statistics[[ref]][["original"]]
  frag.sim=fragment.statistics[[ref]][["simulated"]]
  

  for(enh in c(enhancer.datasets[[ref]], "restriction_fragments")){
    
    enh.obs=enhancer.statistics[[ref]][[enh]][["original"]]
    enh.sim=enhancer.statistics[[ref]][[enh]][["simulated"]]
    
    if (enh == "restriction_fragments"){all=unique(c(frag.obs$ID, frag.sim$ID))}else{all=unique(c(enh.obs$enh, enh.sim$enh))}
    
    path=paste(pathPhyloP, ref, "/", enh, "/phyloP_", way, "_MaskedExons_Ensembl94.txt", sep="")
    
    if(file.exists(path)){
      
      ## we assign a default value to non-lifted elements
      ## default.cons defined in parameters.R
      
      phyloPscore=rep(default.cons, length(all))
      phyloPscore.default0=rep(0, length(all))
      
      names(phyloPscore)=all
      names(phyloPscore.default0)=all

      cons=fread(path, h=T)
      class(cons)="data.frame"
      
      ## use only filtered (non-exonic) sequence covered but at least 10pb
      this.phyloPscore=cons$Score
      this.phyloPscore[which(cons$AnalyzedLength < minAlnLength)] = NA
      names(this.phyloPscore)=cons[,"ID"]
      
      phyloPscore[intersect(names(phyloPscore), names(this.phyloPscore))]=this.phyloPscore[intersect(names(phyloPscore), names(this.phyloPscore))]

      ## default score = 0 for missing base pair
      cons$Score.default0 = cons$Score * (cons$CoveredLength/cons$AnalyzedLength)
      
      this.phyloPscore.0 = cons$Score.default0
      this.phyloPscore.0[which(cons$CoveredLength == 0)] = 0
      this.phyloPscore.0[which(cons$AnalyzedLength < minAlnLength)] = NA
      
      names(this.phyloPscore.0)=cons[,"ID"]
      
      phyloPscore.default0[intersect(names(phyloPscore.default0), names(this.phyloPscore.0))]=this.phyloPscore.0[intersect(names(phyloPscore.default0), names(this.phyloPscore.0))]
      
      ## save results
      save(phyloPscore, phyloPscore.default0, file=paste(pathFigures, "RData/data.phyloP.scores.", enh,".",ref,".RData", sep=""))
      
    }else{print(paste("cannot find file for", ref, enh, sep=" "))}
  }
}

###########################################################################
