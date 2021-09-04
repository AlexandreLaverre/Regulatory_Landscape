############################################################################

## all scripts that generate .RData objects for the figures
## order is important
data.scripts=c()

## 0 dependencies
data.scripts=c(data.scripts, "data.gene.annotations.R")
data.scripts=c(data.scripts, "data.ortho.genes.R")
data.scripts=c(data.scripts, "data.sample.info.R")
data.scripts=c(data.scripts, "data.bait.annotation.R")
data.scripts=c(data.scripts, "data.fragment.statistics.R") ## depends on aberrant fragments but data is re-read
data.scripts=c(data.scripts, "data.Shh.figure.R")
data.scripts=c(data.scripts, "data.sequence.conservation.enhancers.R") ## data is provided for all enhancers 
data.scripts=c(data.scripts, "data.gene.expression.R")

## 1 dependency
data.scripts=c(data.scripts, "data.aberrant.fragments.R") ## depends on sample info
data.scripts=c(data.scripts, "data.enhancer.statistics.R") ## depends on fragment statistics
data.scripts=c(data.scripts, "data.fragment.contacts.R") ## depends on fragment statistics
data.scripts=c(data.scripts, "data.sequence.conservation.fragments.R") ## depends on fragment statistics
data.scripts=c(data.scripts, "data.promoter.enhancer.correlation.R")  ## depends on fragment contacts ## BOOTBCA
data.scripts=c(data.scripts, "data.CM2019.SomaticOrgans.expression.divergence.R") ## depends on gene annotations
data.scripts=c(data.scripts, "data.CM2019.AllOrgans.expression.divergence.R") ## depends on gene annotations

## 2 dependencies
data.scripts=c(data.scripts, "data.sample.clustering.R") ## sample info, fragment contacts
data.scripts=c(data.scripts, "data.enhancer.coverage.R") ## fragment statistics, fragment contacts ## BOOTBCA
data.scripts=c(data.scripts, "data.phyloP.scores.R") ## fragment statistics, enhancer statistics
data.scripts=c(data.scripts, "data.AFC.R") ## fragment contacts, sample info

data.scripts=c(data.scripts, "data.bootstrap.nb.cell.types.R") ## fragment contacts, sample info ## BOOTBCA

## 3 dependencies
data.scripts=c(data.scripts, "data.bootstrap.conservation.distance.repeats.R") ## sequence conservation, enhancer statistics, fragment statistics  ## BOOTBCA
data.scripts=c(data.scripts, "data.bootstrap.phyloP.distance.R") ## enhancer statistics, fragment statistics, sequence conservation ## BOOTBCA
 
## 4 dependencies
data.scripts=c(data.scripts, "data.gene.enhancer.contacts.R") ## gene annotations, sample info,  fragment contacts, enhancer statistics
data.scripts=c(data.scripts, "data.contact.conservation.enhancers.R") ## gene-enhancer contacts, ortho genes, bait annotations,  sequence conservation for enhancers

## 2 dependencies
data.scripts=c(data.scripts, "data.contact.conservation.stats.R") ## sample info, contact conservation enhancers

## 7 dependencies
data.scripts=c(data.scripts, "data.sequence.conservation.stats.R") ## fragment contacts, gene enhancer contacts, enhancer statistics, fragment statistics, sequence conservation fragments, sequence conservation enhancers, phyloP scores

## 2 dependencies, depends on gene enhancer contacts
data.scripts=c(data.scripts, "data.synteny.conservation.R") ## gene-enhancer contacts, sequence conservation for enhancers

## 8 dependencies
data.scripts=c(data.scripts, "data.regland.conservation.R") ## pretty much everything

############################################################################

for(file in data.scripts){
  ## cleanup
  objects=ls()
  objects=setdiff(objects, c("data.scripts", "file"))
  rm(list=objects)

  ## run the script
  print(file)
  source(file)
}

############################################################################


