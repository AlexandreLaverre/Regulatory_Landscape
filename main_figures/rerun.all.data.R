############################################################################

## all scripts that generate .RData objects for the figures

data.scripts=c("data.Shh.figure.R", "data.sample.info.R", "data.sample.clustering.R",  "data.bait.annotations.R", "data.enhancer.coverage.R",  "data.fragment.contacts.R",  "data.gene.enhancer.contacts.R", "data.sequence.conservation.R", "data.sequence.conservation.bycelltype.R", "data.synteny.conservation.R", "data.gene.annotations.R", "data.gene.expression.R", "data.promoter.enhancer.correlation.R")

############################################################################

for(file in data.scripts){
  ## cleanup
  objects=ls()
  objects=setdiff(objects, c("data.scripts", "file"))
  rm(list=objects)

  ## run the script
  source(file)
}

############################################################################


