
## all scripts that generate figures

fig.scripts.fast=c( "Figure2.R", "Figure4.R", "Figure5.R", "Figure6.R")

############################################################################

for(file in fig.scripts.fast){
  ## cleanup
  objects=ls()
  objects=setdiff(objects, c("fig.scripts.fast", "file"))
  rm(list=objects)

  ## run the script
  print(file)
  source(file)
}

############################################################################

fig.scripts.slow=c("Figure1.R", "Figure3.R")

print("Do you want to redo figures 1 and 3 (BCa computations)?")

consent=scan(what=character(), nmax=1)

if(consent=="yes"){

  for(file in fig.scripts.slow){
    ## cleanup
    objects=ls()
    objects=setdiff(objects, c("fig.scripts.slow", "file"))
    rm(list=objects)
    
    ## run the script
    print(file)
    source(file)
  }
}

############################################################################

