#####################################################################

for(file in c("data.gene.annotations.R", "data.bait.annotations.R", "data.interactions.per.sample.R", "data.merged.interactions.R", "data.enhancers.R")){
  objects=setdiff(ls(), file)

  rm(list=objects)

  source(file)
}

#####################################################################
