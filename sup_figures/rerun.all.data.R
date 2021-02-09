############################################################################

## all scripts that generate .RData objects for the figures

data.scripts=c("data.features.coverage.R", "data.common.cell.expression.divergence.R", "data.common.cell.regland.conservation.R", "data.cell.types.parallel.trends.R", "data.samples.cumulative.interactions.R")

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


