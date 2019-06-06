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

# viruses = c("H1N1pandemic", "H1N1seasonal", "H3N2", "H2N2", "InfB")
# virus_labels = c("p09 like H1N1", "pre 09 H1N1", "H3N2", "H2N2", "Influenza B")
# virusname = c("h1n1pdm", "h1n1sea", "h3n2", "h2n2", "infB")
# 
viruses = c("H3N2", "H3N2", "H3N2", "H3N2", "H3N2")
virus_labels = c("H3N2 1980-2020", "H3N2 1980-1990", "H3N2 1990-2000", "H3N2 2000-2010", "H3N2 2010-2017")
virusname = c("h3n2", "h3n2ancient", "h3n2old","h3n2recent","h3n2new")


# viruses = c("H1N1seasonal", "H1N1pandemic" , "H3N2", "H2N2", "InfB")
# virus_labels = c("pre 09 H1N1", "p09 like H1N1",  "H3N2", "H2N2", "Influenza B")
# virusname = c("h1n1sea", "h1n1pdm", "h3n2", "h2n2", "infB")


first = T;
for(i in seq(1, length(viruses))){
  # get the names run log files
  # set the log file name
  filename = path=paste("./combined/", virusname[[i]],".trunk.txt", sep="")
  
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
ggsave(plot=p,paste("../../../Reassortment-Text/Figures/trunk_rates_h3n2.pdf", sep=""),width=10, height=3)
