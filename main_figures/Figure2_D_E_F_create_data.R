library(Hmisc)
options(stringsAsFactors = FALSE)

ref_sp = "mouse"

path <- "/home/laverre/Manuscript"

################################## Correlation Gene expression & Nb enhancers ################################## 
## Gene expression
expdiv=read.table(paste(path, "/SupplementaryDataset6/ExpressionDivergence_CardosoMoreira2019.txt", sep=""), h=T, stringsAsFactors=F, sep="\t")

if(ref_sp == "human"){
  rownames(expdiv)=expdiv$IDHuman
  enhancers <- c("CAGE", "ENCODE", "RoadMap", "GRO_seq")

  }else{
    rownames(expdiv)=expdiv$IDMouse
    enhancers <- c("CAGE", "ENCODE")
    }

## Gene Hi-C contacts

if (ref_sp == "mouse"){nb_class_CAGE=12}else{nb_class_CAGE=11}
gene_expression_enhancers <- data.frame(matrix(vector(), 10, 1))

for (enh in enhancers){
  print(enh)
  regland = read.table(paste(path, "/SupplementaryDataset7/", ref_sp, "/evolution_summary_by_gene/", enh, "/", enh, "_original_summary_conserv_all_0.5.txt",sep=""), h=T, sep="\t", row.names = 1)
  regland <- regland[which(regland$nb_total < 250),]
  common=intersect(rownames(expdiv), rownames(regland))
  expdiv=expdiv[common,]
  regland=regland[common,]
  
  # Made quantile of nb enhancers
  if(enh=="CAGE"){regland$fraction_group <- cut2(regland$nb_total, g=nb_class_CAGE)
  }else{regland$fraction_group <- cut2(regland$nb_total, g=10)}
  
  
  box <- boxplot(log2(expdiv$Human_MeanRPKM+1)~regland$fraction_group, plot=F)
  gene_expression_enhancers[[paste0(enh)]] <- box$stats[3,]
  gene_expression_enhancers[[paste0(enh, "_conflow")]] <-  box$conf[1,]
  gene_expression_enhancers[[paste0(enh, "_confup")]] <- box$conf[2,]
  
  R=cor.test(regland$nb_total,log2(expdiv$MeanRPKM+1), method="pearson")
  print(R)
  rho=cor.test(regland$nb_total,log2(expdiv$MeanRPKM+1), method="spearman")
  print(rho)
  
}

################################## Correlation Gene Expression & Enhancers activity ################################## 
path_correl <- paste(path, "SupplementaryDataset8", sep="/")
minDistance=25e3
maxDistance=2.5e6


for (enh in enhancers){
  obs <- read.table(paste(path_correl, ref_sp, enh, "expression_correlations_real_data.txt", sep="/"), h=T, sep="\t")
  simul <- read.table(paste(path_correl, ref_sp, enh, "expression_correlations_simulated_data.txt", sep="/"), h=T, sep="\t")

  # According to distance
  obs$dist_class <-cut(obs$Distance, breaks=seq(from=minDistance, to=maxDistance+50000, by=50000), include.lowest = T)
  simul$dist_class <- cut(simul$Distance, breaks=seq(from=minDistance, to=maxDistance+50000, by=50000), include.lowest = T)
  
  if(enh == "CAGE"){
    obs_correl_activity_dist <- data.frame(matrix(vector(), length(levels(obs$dist_class)), 1))
    simul_correl_activity_dist <- data.frame(matrix(vector(), length(levels(simul$dist_class)), 1))
  }
  
  obs_correl_activity_dist[[paste0(enh)]] <- sapply(levels(obs$dist_class), function(x)
    mean(obs[which(obs$dist_class == x),]$SpearmanCorrelation, na.rm=T))
  
  obs_correl_activity_dist[[paste0(enh, "_conflow")]] <- sapply(levels(obs$dist_class), function(x)  
    t.test(obs[which(obs$dist_class == x),]$SpearmanCorrelation, na.rm=T)[["conf.int"]][1])
  
  obs_correl_activity_dist[[paste0(enh, "_confup")]] <- sapply(levels(obs$dist_class), function(x)  
    t.test(obs[which(obs$dist_class == x),]$SpearmanCorrelation, na.rm=T)[["conf.int"]][2])
  
  
  simul_correl_activity_dist[[paste0(enh)]] <- sapply(levels(simul$dist_class), function(x)
    mean(simul[which(simul$dist_class == x),]$SpearmanCorrelation))
  
  simul_correl_activity_dist[[paste0(enh, "_conflow")]] <- sapply(levels(simul$dist_class), function(x)  
    t.test(simul[which(simul$dist_class == x),]$SpearmanCorrelation)[["conf.int"]][1])
  
  simul_correl_activity_dist[[paste0(enh, "_confup")]] <- sapply(levels(simul$dist_class), function(x)  
    t.test(simul[which(simul$dist_class == x),]$SpearmanCorrelation)[["conf.int"]][2])
  

    
  }
  
correl_activity <- list(obs=obs_correl_activity_dist, simul=simul_correl_activity_dist)

save(gene_expression_enhancers, correl_activity, file = paste(path, "/Figures/Fig2_", ref_sp, "_D_E.Rdata", sep=""))

  
