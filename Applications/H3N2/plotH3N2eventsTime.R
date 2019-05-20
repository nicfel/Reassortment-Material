######################################################
######################################################
# Here the inferred mean coalescent and migration
# rate ratios are plotted
######################################################
######################################################
library(ggplot2)
library(colorblindr)

# clear workspace
rm(list = ls())

# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

segments = c("HA", "MP", "NA", "NP", "NS1", "PA", "PB1", "PB2", "prior")

threshold = 0.5

first = T;
for (a in seq(1,length(segments)-2)){
  for (b in seq(a+1,length(segments)-1)){
    # get the names run log files
    filename = paste(path="./combined/h3n2.event.",segments[[a]], "_", segments[[b]], ".txt",sep="")
    
    t <- read.table(filename, header=TRUE, sep="\t")
    t$segment.left = segments[[a]]
    t$segment.right = segments[[b]]

    if (first){
      events=t
      first=F
    }else{
      events=rbind(events,t)
    }

  }
}

# rate$seg1 <- factor(rate$seg1, levels = c("HA", "MP", "NA", "NP", "NS1", "PA", "PB1", "PB2"))
# rate$seg2 <- factor(rate$seg2, levels = c("HA", "MP", "NA", "NP", "NS1", "PA", "PB1", "PB2"))

p.rea.count <- ggplot(events) +
  geom_density(aes(2017-height.mean, weight=posterior)) +
  facet_grid(segment.left~segment.right)
  xlab("") +
  ylab("") +
  theme_light() 

plot(p.rea.count)

ggsave(plot=p.rea.count,paste("../../../Reassortment-Text/Figures/pairwise/h1n1_count", ".pdf" , sep=""),width=10, height=3)


