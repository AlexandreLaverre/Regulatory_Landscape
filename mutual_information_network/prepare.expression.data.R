#########################################################################################

path="/sps/biometr/necsulea/RegulatoryLandscapes/"
pathData=paste(path, "data/FANTOM5/", sep="")
pathResults=paste(path, "results/mutual_information_network/", sep="")

#########################################################################################

genomes=c("hg19", "mm9")
species=c("human", "mouse")

minvalue=1
minsamples=10

#########################################################################################

for(i in 1:2){
  sp=species[i]
  genome=genomes[i]

  print("reading gene data")
  exp.genes=read.table(paste(pathData, genome, "/", genome, ".gene_phase1and2combined_tpm.osc.averagevalues.txt", sep=""), h=T, stringsAsFactors=F, sep="\t")
  print("done")

  print("reading enhancer data")
  exp.enhancers=read.table(paste(pathData, genome, "/", sp, "_permissive_enhancers_phase_1_and_2_expression_tpm_matrix.averagevalues.txt", sep=""), h=T, stringsAsFactors=F, sep="\t")
  print("done")

  print(all(colnames(exp.genes)==colnames(exp.enhancers)))

  samples=colnames(exp.genes)[-1]

  nbsamples.genes=apply(exp.genes[,samples],1, function(x) length(which(x>=minvalue)))
  nbsamples.enhancers=apply(exp.enhancers[,samples],1, function(x) length(which(x>=minvalue)))

  exp.selected.genes=exp.genes[which(nbsamples.genes>=minsamples),]
  exp.selected.enhancers=exp.enhancers[which(nbsamples.enhancers>=minsamples),]

  writeLines(exp.selected.enhancers$ID, con=paste(pathResults, genome, "/EnhancerList.txt", sep=""))

  exp.data=rbind(exp.selected.genes, exp.selected.enhancers, stringsAsFactors=F)
  colnames(exp.data)[1]="Gene"

  write.table(exp.data, file=paste(pathResults, genome, "/TPM.txt", sep=""), row.names=F, col.names=T, sep="\t", quote=F)
  
}

#########################################################################################
