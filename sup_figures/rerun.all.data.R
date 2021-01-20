############################################################################

## all scripts that generate .RData objects for the figures

all.files=system("ls", intern=T)
data.scripts=grep("^data", all.files, value=T)

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


