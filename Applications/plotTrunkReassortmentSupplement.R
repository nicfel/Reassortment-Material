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

viruses = c("H3N2", "H3N2", "H3N2", "InfB", "InfB", "InfB")
virus_labels = c("H3N2 2y", "H3N2 4y", "H3N2 6y", "InfB 2y", "InfB 4y", "InfB 6y")
virusname = c("h3n2.2", "h3n2.", "h3n2.6", "infB.2", "infB.", "infB.6")



first = T;
for(i in seq(1, length(viruses))){
  # get the names run log files
  # set the log file name
  filename = path=paste("./",viruses[[i]] ,"/combined/", virusname[[i]],"trunk.txt", sep="")
  
  # read the file in
  t <- read.table(filename, header=F, sep="\t")
  
  new.rate.diff = data.frame(virus=virus_labels[[i]], rate.diff=t$V1/t$V3-t$V2/t$V4)
  
  new.trunk.rate = data.frame(rate=t$V1/t$V3, virus = paste(virus_labels[[i]]), trunk="on trunk rates", ratio="number of events per lineage and year")
  new.trunk.rate = rbind(new.trunk.rate, data.frame(rate=t$V2/t$V4,  virus = paste(virus_labels[[i]]), trunk="off trunk rates", ratio="number of events per lineage and year"))
  new.trunk.rate = rbind(new.trunk.rate, data.frame(rate=t$V1/t$V3-t$V2/t$V4,  virus = paste(virus_labels[[i]]), trunk="rate difference", ratio="difference in the number of events per lineage and year"))
  
  if (first){
    rate.diff = new.rate.diff
    trunk.rate = new.trunk.rate
    first = F;
  }else{
    rate.diff = rbind(rate.diff, new.rate.diff)
    trunk.rate = rbind(trunk.rate, new.trunk.rate)
  }
  
}

# segments = c("HA", "MP", "NA", "NP", "NS1", "PA", "PB1", "PB2")

p <- ggplot(rate.diff) +
  geom_violin(aes(x=virus,y=rate.diff))+
  theme_minimal()

plot(p)

p <- ggplot(trunk.rate) +
  geom_violin(aes(x=virus,y=rate, fill=trunk))+
  scale_fill_OkabeIto() +
  ylab("") +
  theme(legend.position = "none")+
  facet_wrap(ratio~., scale="free_y", ncol=1) +
  theme_minimal() +
  theme(strip.background = element_blank(),
  strip.text.y = element_blank()) +
  xlab("")

plot(p)
ggsave(plot=p,paste("../../Reassortment-Text/Figures/trunk_rates_distance.pdf", sep=""),width=10, height=3)
