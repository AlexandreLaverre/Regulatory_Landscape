############################################################################

## all scripts that generate .RData objects for the figures
## we need a certain order
data.scripts=c("data.sample.info.R", "data.aberrant.fragments.R", "data.fragment.statistics.R", "data.enhancer.statistics.R", "data.bait.annotation.R",  "data.fragment.contacts.R", "data.gene.enhancer.contacts.R",  "data.Shh.figure.R",  "data.sample.clustering.R", "data.enhancer.coverage.R", "data.sequence.conservation.R", "data.synteny.conservation.R", "data.gene.annotations.R", "data.gene.expression.R", "data.promoter.enhancer.correlation.R", "data.synteny.conservation.R", "data.contact.conservation.enhancers.R", "data.contact.conservation.stats.R", "data.AFC.R", "data.common.cell.regland.conservation.R", "data.common.cell.expression.divergence.R",  "data.regland.conservation.R", "data.CM2019.SomaticOrgans.expression.divergence.R")

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


