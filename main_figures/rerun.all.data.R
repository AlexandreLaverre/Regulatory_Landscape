############################################################################

## all scripts that generate .RData objects for the figures

data.scripts=c("data.fragment.statistics.R", "data.enhancer.statistics.R", "data.bait.annotation.R", "data.gene.enhancer.contacts.R", "data.fragment.contacts.R",  "data.Shh.figure.R", "data.sample.info.R", "data.sample.clustering.R",  "data.enhancer.coverage.R",   "data.sequence.conservation.R", "data.sequence.conservation.bycelltype.R", "data.synteny.conservation.R", "data.gene.annotations.R", "data.gene.expression.R", "data.promoter.enhancer.correlation.R", "data.contact.conservation.enhancers.R", "data.contact.conservation.enhancers.stats.R", "data.AFC.R")

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


