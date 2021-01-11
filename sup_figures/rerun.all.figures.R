
## all scripts that generate figures

fig.scripts=system("ls", intern=T)
fig.scripts=grep("SupplementaryFigure", fig.scripts, value=T)
fig.scripts=grep("SupplementaryFigureX", fig.scripts, value=T, invert=T)

############################################################################

for(file in fig.scripts){
  ## cleanup
  objects=ls()
  objects=setdiff(objects, c("fig.scripts", "file"))
  rm(list=objects)

  ## run the script
  print(file)
  source(file)
}

############################################################################


