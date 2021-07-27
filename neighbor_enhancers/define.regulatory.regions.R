##########################################################################

pathResults="../../results/neighbor_enhancers/"

release=94

##########################################################################

for(sp in c("human", "mouse")){
  tss=read.table(paste(pathResults, tolower(sp), "/canonical_transcripts_Ensembl",release,".txt",sep=""), h=T, stringsAsFactors=F, sep="\t")

  regions=c()

  for(chr in c(as.character(1:22), "X", "Y")){
    print(chr)
    
    this.tss=tss[which(tss$chr==chr),]

    if(nrow(this.tss)>0){
      
      this.tss=this.tss[order(this.tss$TSS),]
      
      ## we initialize regulatory region - 5kb upstream, 1kb downstream
      
      this.tss$start_region=rep(NA, dim(this.tss)[1])
      this.tss$end_region=rep(NA, dim(this.tss)[1])
      
      fwd=which(this.tss$strand==1)
      rev=which(this.tss$strand==-1)
      
      this.tss$start_region[fwd]=this.tss$TSS[fwd]-5000
      this.tss$end_region[fwd]=this.tss$TSS[fwd]+1000
      
      this.tss$start_region[rev]=this.tss$TSS[rev]-1000
      this.tss$end_region[rev]=this.tss$TSS[rev]+5000
      
      for(i in 1:nrow(this.tss)){
        if(i>1){
          previousTSS=this.tss$TSS[i-1]
          this.tss$start_region[i]=min(this.tss$start_region[i],  max(this.tss$TSS[i]-1e6, previousTSS+1))
        } else{
          this.tss$start_region[i]=max(this.tss$TSS[i]-1e6, 1)
        }
        
        if(i<nrow(this.tss)){
          nextTSS=this.tss$TSS[i+1]
          this.tss$end_region[i]=max(this.tss$end_region[i], min(this.tss$TSS[i]+1e6, nextTSS-1))
        } else{
          this.tss$end_region[i]=this.tss$TSS[i]+1e6
        }
      }
      
      if(length(regions)==0){
        regions=this.tss
      } else{
        regions=rbind(regions, this.tss)
      }
    }
  }
  
  write.table(regions, file=paste(pathResults, tolower(sp), "/regulatory_regions_Ensembl",release,".txt", sep=""), row.names=F, col.names=T, sep="\t", quote=F)
}

##########################################################################
