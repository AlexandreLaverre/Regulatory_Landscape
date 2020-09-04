#######################################################################

source("parameters.R")

load(paste(pathFigures, "RData/data.contact.conservation.RData", sep=""))
load(paste(pathFigures, "RData/data.sample.info.RData", sep=""))

#######################################################################

for(ref in c("human", "mouse")){
  tg=setdiff(c("human", "mouse"), ref)

  ref.samples=toupper(sampleinfo[[ref]][,"Acronym"])
  tg.samples=toupper(sampleinfo[[tg]][,"Acronym"])
  
  for(enh in enhancer.datasets[[ref]]){
    cons.real=contact.conservation[[paste(ref, "2", tg, sep="")]][[enh]][["obsobs"]]
    cons.sim=contact.conservation[[paste(ref, "2", tg, sep="")]][[enh]][["obssim"]]

    cons.real$ID=paste(cons.real$origin_gene, cons.real$origin_enh, sep="_")
    cons.sim$ID=paste(cons.sim$origin_gene, cons.sim$origin_enh, sep="_")

    colnames(cons.real)=toupper(colnames(cons.real))
    colnames(cons.sim)=toupper(colnames(cons.sim))

    consmat.real=matrix(rep(NA, length(ref.samples)*length(tg.samples)), nrow=length(ref.samples))
    rownames(consmat.real)=ref.samples
    colnames(consmat.real)=tg.samples

    consmat.sim=consmat.real

    for(rs in ref.samples){
      this.real=cons.real[which(cons.real[,rs]!=0),]
      this.sim=cons.sim[which(cons.sim[,rs]!=0),]
      
      for(ts in tg.samples){
        fr.real=length(unique(this.real$ID[which(this.real[,ts]!=0)]))/length(unique(this.real$ID))
        consmat.real[rs,ts]=fr.real

        fr.sim=length(unique(this.sim$ID[which(this.sim[,ts]!=0)]))/length(unique(this.sim$ID))
        consmat.sim[rs,ts]=fr.sim
       
      }
    }
  }

}

#######################################################################

## duplications dans le tableau de contact!
## probablement cas ou enhancer a cheval sur plusieurs fragments
## comment calcule-t-on la conservation pour la simulation?
## donnees reelles ref vs. donnees simulees cible 

#######################################################################
