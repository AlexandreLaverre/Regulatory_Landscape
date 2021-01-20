
## all scripts that generate figures

list=system("ls", intern=T)
fig.scripts=grep("SupplementaryFigure", list, value=T)
fig.scripts=grep("SupplementaryFigureX", fig.scripts, value=T, invert=T)
extended.fig.scripts=grep("Extended", list, value=T)

fig.scripts = append(extended.fig.scripts, fig.scripts)
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


