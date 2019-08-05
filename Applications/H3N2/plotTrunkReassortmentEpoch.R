######################################################
######################################################
# Here the inferred mean coalescent and migration
# rate ratios are plotted
######################################################
######################################################
library(ggplot2)
library(colorblindr)
library(grid)
library(gridExtra)

# clear workspace
rm(list = ls())

# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

viruses = c("H3N2", "H3N2", "H3N2","H3N2", "H3N2", "H3N2","H3N2", "H3N2", "H3N2")

virus_labels = c("H3N2 1980-2019","H3N2 1980-2019","H3N2 1980-2019",
                 "H3N2 1980-2005","H3N2 1980-2005","H3N2 1980-2005",
                 "H3N2 2005-2019","H3N2 2005-2019","H3N2 2005-2019")

virusname = c("h3n2_sub1","h3n2_sub1","h3n2_sub1", 
              "h3n2old","h3n2old","h3n2old",
              "h3n2recent","h3n2recent","h3n2recent")

trunk = c(2,4,6,
          2,4,6,
          2,4,6)

replicate = c("1","2","3",
              "1","2","3", 
              "1", "2", "3", 
              "2", 
              "1", "2", "3")

offset=0.1
range = 0.35

xval = c(1*range-offset,1*range,1*range+offset,
         2*range-offset,2*range,2*range+offset,
         3*range-offset,3*range,3*range+offset)
         

first = T;
for(i in seq(1, length(viruses))){
  # get the names run log files
  # set the log file name
  filename = path=paste("./combined/", virusname[[i]],".", trunk[[i]],"trunk.txt", sep="")
  
  # read the file in
  t <- read.table(filename, header=F, sep="\t")
  
  new.rate.diff = data.frame(virus=virus_labels[[i]], x=xval[i], rate.diff=t$V1/t$V3-t$V2/t$V4)
  
  new.trunk.rate = data.frame(rate=t$V1/t$V3, x=xval[i], virus = paste(virus_labels[[i]]), trunk="fit", ratio="absolute values", replicate=replicate[[i]])
  new.trunk.rate = rbind(new.trunk.rate, data.frame(rate=t$V2/t$V4,x=xval[i],   virus = paste(virus_labels[[i]]), trunk="unfit", ratio="absolute values", replicate=replicate[[i]]))
  
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

p.fit <- ggplot(trunk.rate) +
  geom_violin(aes(x=x,y=rate,group=interaction(x, trunk), fill=trunk, color=trunk))+
  scale_fill_OkabeIto(name="") +
  scale_color_OkabeIto(name="") +
  ylab("reassortment rate\nper lineage and year") +
  theme(legend.position = "none")+
  theme_minimal() +
  theme(strip.background = element_blank(),
        legend.position = c(0.9,0.8),
        strip.text.y = element_blank()) +
  scale_x_continuous(breaks=seq(1,3)*range, labels=c("H3N2 1980-2019", "H3N2 1980-2005", "H3N2 2005-2019"))+
  xlab("")

p.diff <- ggplot(rate.diff) +
  geom_violin(aes(x=x,y=rate.diff,group=x))+
  scale_fill_OkabeIto(name="") +
  scale_color_OkabeIto(name="") +
  ylab("reassortment rate difference\nper lineage and year") +
  theme(legend.position = "none")+
  theme_minimal() +
  theme(strip.background = element_blank(),
        strip.text.y = element_blank()) +
  scale_x_continuous(breaks=seq(1,3)*range, labels=c("H3N2 1980-2019", "H3N2 1980-2005", "H3N2 2005-2019"))+
  xlab("") +
  geom_hline(yintercept=0)

plot.all = grid.arrange(p.fit, p.diff,ncol=1)


ggsave(plot=plot.all,paste("../../../Reassortment-Text/Figures/trunk_rates_h3n2_epoch.pdf", sep=""),width=7, height=5)
