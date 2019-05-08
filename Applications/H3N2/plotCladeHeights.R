######################################################
######################################################
# combine the gene tree runs and run the mcc trees
######################################################
######################################################
# clear workspace
library(ggplot2)
library(grid)
library(gridExtra)
rm(list = ls())

# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

heights <- list.files(path=paste("./trees/heights",sep=""), pattern="*\\.txt$", full.names = TRUE)

t = read.table("./trees/Heights.csv", header=T, sep=",", na.strings="")

t$mid1 = (t$lower1+t$upper1)/2
t$mid2 = (t$lower2+t$upper2)/2
t$ratio = (t$upper1-t$lower1)/(t$upper2-t$lower2)
t = t[-which(t$post1<0.25 | t$post2<0.25),]


unique.times = unique(t$time)
unique.segment = unique(t$segment)
p <- list()
c = 1
for (a in seq(1,length(unique.times))){
  for (b in seq(1,length(unique.segment))){
    subset = t[which(t$segment==unique.segment[[b]] & t$time==unique.times[[a]]),]
    p[[c]] <- ggplot(data=subset) +
      geom_histogram(aes(ratio)) +
      geom_vline(xintercept=mean(subset$ratio), color="red")+
      scale_x_log10(limits=c(0.1,10)) +
      ggtitle(paste(unique.segment[[b]], unique.times[[a]], "mean =", mean(subset$ratio)))
    
    c=c+1
  }
}

plot.sampling = do.call("grid.arrange",c(p, ncol=8))
plot(plot.sampling)
