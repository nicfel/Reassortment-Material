######################################################
######################################################
# Here the inferred mean coalescent and migration
# rate ratios are plotted
######################################################
######################################################
library(ggplot2)
library(coda)
library(colorblindr)

# clear workspace
rm(list = ls())

# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

# define the threshold of when to count an event
threshold = 0

# define the number of segments
nrSegments = 9

# read in the file with the reassortment distances
fileName="combined/h1n1sea.distance.txt"
con=file(fileName,open="r")
line=readLines(con) 
first = T
for (i in seq(1,length(line))){
  print(paste(i, "of", length(line)) )
  
  splitline = strsplit(line[[i]], split="\t")
  firstevent=T
  remove(dat)
  
  for (j in seq(2,length(splitline[[1]]))){
    tmp = strsplit(splitline[[1]][[j]], split=",")
    for (k in seq(2,length(tmp[[1]]))){
      tmp2 = strsplit(tmp[[1]][[k]], split=":")
    }
    
    tmp2 = strsplit(tmp[[1]][[1]], split=":")
    val1 = gsub("\\{", "",tmp2[[1]][[2]])
    val1 = gsub("\\}", "",val1)
    val2 = gsub("\\{", "",tmp2[[1]][[3]])
    val2 = gsub("\\}", "",val2)
    
    val1 = gsub("8", "",val1)
    val2 = gsub("8", "",val2)
    
    l1 = length(strsplit(val1, split="  ")[[1]])
    l2 = length(strsplit(val2, split="  ")[[1]])
    
    new.new.event = data.frame(isTrunk=as.logical(tmp2[[1]][[1]]), seg1=max(l1,l2), seg2 = min(l1,l2))
    if (firstevent){
      new.event = new.new.event
      firstevent = F
    }else{
      new.event = rbind(new.event, new.new.event)
    }
  }
  
  indices_trunk = which(new.event$isTrunk)
  indices_off = which(!new.event$isTrunk)
  
  mean1 = mean(new.event[indices_trunk, "seg2"])
  mean2 = mean(new.event[indices_off, "seg2"])
  
  if (mean1<mean2){
    new.nr.event = data.frame(larger=1)
  }else{
    new.nr.event = data.frame(larger=0)
  }
  # get the distribution of the number of events
  if (first){
    nr.event = new.nr.event
    first = F
  }else{
    nr.event = rbind(nr.event, new.nr.event)
  }
}
