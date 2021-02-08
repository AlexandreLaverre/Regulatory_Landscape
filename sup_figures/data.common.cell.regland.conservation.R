#######################################################################################

options(stringsAsFactors = FALSE)

source("../main_figures/parameters.R")

#######################################################################################

print("loading data")

load(paste(pathFigures,"RData/data.sample.info.RData", sep=""))
load(paste(pathFigures, "RData/data.gene.enhancer.contacts.RData", sep=""))
load(paste(pathFigures, "RData/data.ortho.genes.RData", sep=""))
load(paste(pathFigures, "RData/data.contact.conservation.enhancers.RData", sep=""))
load(paste(pathFigures, "RData/data.bait.annotation.RData", sep=""))
load(paste(pathFigures, "RData/data.fragment.contacts.RData",sep=""))

print("done")

## we load sequence and synteny conservation later

#######################################################################################

## min and max number of enhancers for contacts and synteny conservation

min.nb=5
max.nb=100

#######################################################################################

## extract samples for common cell types

cells <- c("embryonic stem cells", "pre-adipocytes", "B lymphocytes")
syncells <- c("ESC", "adipo", "Bcell")

common.cells <- lapply(cells, function(x) lapply(c("human", "mouse"), function(y) {z=sampleinfo[[y]]; z[which(z$Broad.cell.type.or.tissue==x),"Sample.ID"]}))

names(common.cells)=syncells

for(c in syncells){
  names(common.cells[[c]])=c("human", "mouse")
}

#######################################################################################

regland.conservation=list()

for(ref in c("human", "mouse")){
  print(ref)
  ## prepare results
  regland.conservation[[ref]]=list()
  
  tg=setdiff(c("human", "mouse"), ref)
  
  for(enh in enhancer.datasets[[ref]]){
    print(enh)

    ## prepare results
    regland.conservation[[ref]][[enh]]=list()
    
    ## we start from the gene-enhancer contact data
    
    contacts=gene.enhancer.contacts[[ref]][[enh]][["real"]]
    
    ## we select only previously filtered orthologous genes
    
    contacts=contacts[which(contacts$gene%in%ortho[,ref]),]
    
    ## we add sequence conservation column

    load(paste(pathFigures, "RData/data.sequence.conservation.enhancers.",enh,".",ref,"2", tg,".RData", sep=""))
    contacts$align_score=pcungapped[contacts$enhancer]

    ## we load synteny conservation - already filtered, min sequence conservation >= 10% threshold
    
    load(paste(pathFigures, "/RData/data.synteny.conservation.", ref, ".RData", sep=""))
    
    synteny=conserv_synteny[[enh]][[tg]][["synt_obs"]]
    
    ## we load contact conservation - already filtered, min sequence conservation 0.4

    contact.cons=contact.conservation[[paste(ref, "2", tg, sep="")]][[enh]][["obs"]]
    
    for (cell in syncells){
      print(cell)
            
      ## Take only contacts present in reference cell
      filtered.contacts <- contacts[which(apply(contacts[common.cells[[cell]][[ref]]], 1, function(x) any(x>0))),]

      ## prepare results      
      all.genes = levels(factor(filtered.contacts$gene))
      results=list("gene"=all.genes)
            
      ## nb enhancers by gene

      nb.contacts=as.numeric(table(factor(filtered.contacts$gene, levels=all.genes)))

      ## divide nb of contacts into 5 classes 
     
      results[["nb.contacts"]]=nb.contacts
      results[["class.nb.contacts"]]=cut(nb.contacts, breaks=c(0, 5, 10, 15, 20, max(nb.contacts)), include.lowest=T, labels=c("1-5", "6-10", "11-15", "16-20", ">20"))
      
      ## mean conservation score by gene
      
      mean.aln.score=tapply(filtered.contacts$align_score, factor(filtered.contacts$gene, levels=all.genes), mean, na.rm=T)
      
      results[["mean.aln.score"]]=mean.aln.score
      results[["class.aln.score"]]=cut(mean.aln.score, breaks=seq(from=0, to=1, length=6), include.lowest=T)
      
      ## synteny conservation
      
      filtered.synteny=synteny[which(synteny$id%in%filtered.contacts$id),]
      
      fr.cons.synt=tapply(filtered.synteny$target_dist,  factor(filtered.synteny$origin_gene, levels=all.genes), function(x) length(which(x<maxDistanceSyntenyTarget))/length(x))
      
      nb.enh.synt=as.numeric(table(factor(filtered.synteny$origin_gene, levels=all.genes)))
      
      ## if fewer than min.nb enhancers, we assign NA values to synteny conservation
      
      fr.cons.synt[which(nb.enh.synt < min.nb | nb.enh.synt > max.nb)]=NA
      
      results[["fr.synteny.cons"]]=fr.cons.synt
      
      results[["class.synteny.cons"]]=factor(fr.cons.synt<1,levels=c(TRUE, FALSE), labels=c("broken", "conserved"))
      
      ## contact conservation
      
      filtered.contact.cons=contact.cons[which(contact.cons$id%in%filtered.contacts$id),]
      
      filtered.contact.cons$is_conserved=apply(filtered.contact.cons[common.cells[[cell]][[tg]]], 1, function(x) any(x>0))
      
      fr.cons.contact=tapply(filtered.contact.cons$is_conserved, factor(filtered.contact.cons$origin_gene, levels=all.genes), function(x) length(which(x))/length(x))
      
      nb.enh.contact=as.numeric(table(factor(filtered.contact.cons$origin_gene, levels=all.genes)))
      
      ## if fewer than min.nb enhancers, we assign NA values to synteny conservation
      
      fr.cons.contact[which(nb.enh.contact < min.nb | nb.enh.contact > max.nb)]=NA
      
      results[["fr.contact.cons"]]=fr.cons.contact
      
      results[["class.contact.cons"]]=cut(fr.cons.contact, breaks=c(0, 0.1, 0.4, 1), include.lowest=T)
      
      ## format results
      results=as.data.frame(results)
      
      regland.conservation[[ref]][[enh]][[cell]]=results
    }
  }
}

####################################################################################################################

## save results

save(regland.conservation, file=paste(pathFigures, "RData/data.common.cells.regland.conservation.RData",sep=""))

####################################################################################################################
