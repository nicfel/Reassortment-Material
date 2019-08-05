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

first = T

# read in the file with the reassortment distances
fileName= paste("./combined/h3n2new.distance.txt",sep="")
con=file(fileName,open="r")
line=readLines(con) 
close(con)
for (i in seq(1,length(line))){
  print(paste(i, "of", length(line)) )
  
  splitline = strsplit(line[[i]], split="\t")
  firstevent=T
  remove(dat)
  c=1
  
  for (j in seq(2,length(splitline[[1]]))){
    tmp = strsplit(splitline[[1]][[j]], split=",")
    tmp2 = strsplit(tmp[[1]][[1]], split=":")
    # checks if the event is on the trunk
    isTrunk = as.logical(tmp2[[1]][[1]])
    l1 = length(strsplit(tmp2[[1]][[2]], split="  ")[[1]])
    l2 = length(strsplit(tmp2[[1]][[3]], split="  ")[[1]])
    randval = rbinom(1,(l1+l2), 0.5 )
    while (randval==0 || randval==(l1+l2)){
      randval = rbinom(1,(l1+l2), 0.5 )
    }
    new.data = data.frame(seg1=min(l1,l2)/(l1+l2), random=min(randval,(l1+l2-randval))/(l1+l2), trunk = isTrunk)
    
    if (firstevent){
      firstevent=F
      data = new.data
    }else{
      data = rbind(data, new.data)
    }
  }
  
  new.reascount = data.frame(larger=1, trunk=mean(data[which(data$trunk),"seg1"]), nontrunk=mean(data[which(!data$trunk),"seg1"]),
                             prior.trunk=mean(data[which(data$trunk),"random"]), prior.nontrunk=mean(data[which(!data$trunk),"random"]))
  
  if(first){
    reascount = new.reascount
    first = F
  }else{
    reascount = rbind(reascount, new.reascount)
  }
  
}


p <- ggplot(reascount) +
  geom_violin(aes(x="ratio",y=trunk/prior.trunk)) +
  geom_violin(aes(x="prior ratio",y=nontrunk/prior.nontrunk)) +
  theme_light() +
  scale_fill_OkabeIto()
  
plot(p)
