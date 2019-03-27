######################################################
######################################################
# Here the inferred mean coalescent and migration
# rate ratios are plotted
######################################################
######################################################
library(ggplot2)

# clear workspace
rm(list = ls())

# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)


# use the matlab standard colors to plot
col0 <- rgb(red=0.0, green=0.4470,blue=0.7410)
col1 <- rgb(red=0.8500, green=0.3250,blue=0.0980)
col2 <- rgb(red=0.9290, green=0.6940,blue=0.1250)
col4 <- rgb(red=0.4660, green=0.6740,blue=0.1880)
col3 <- rgb(red=0.3010, green=0.7450,blue=0.9330)

viruses = c("H1N1pandemic", "H1N1seasonal", "H3N2", "H5N1", "InfB", "InfC", "InfD")

virus_labels = c("p2009 like H1N1", "seasonal H1N1", "H3N2", "H5N1", "Influenza B yam + vic", "Influenza C", "Influenza D")


first = T;
for(i in seq(1, length(viruses))){
  # get the names run log files
  log <- list.files(path=paste("./",viruses[[i]] ,"/out", sep=""), pattern="*[012].log", full.names = TRUE)
  # read them in
  for (j in seq(1, length(log))){
    t <- read.table(log[j], header=TRUE, sep="\t")
    # take a 10 % burn in
    t <- t[-seq(1,ceiling(length(t$Sample)/2)), ]
    
    new.reassortment = data.frame(rate=t$reassortmentRate,Neff=t$popSize.t)
    new.reassortment$virus = virus_labels[[i]]
    if (first){
      reassortment = new.reassortment
      first = F;
    }else{
      reassortment = rbind(reassortment, new.reassortment)
    }
  }
}

p.rea <- ggplot(reassortment) +
  geom_violin(aes(x=virus, y=rate)) +
  scale_y_log10() +
  xlab("") +
  ylab("reassortment rate") +
  theme_minimal()

p.ne <- ggplot(reassortment) +
  geom_violin(aes(x=virus, y=Neff)) +
  scale_y_log10() +
  xlab("") +
  ylab("effective population size") +
  theme_minimal()

p.ratio <- ggplot(reassortment) +
  geom_violin(aes(x=virus, y=rate/Neff)) +
  scale_y_log10() +
  xlab("") +
  ylab("reassortment rate over Ne") +
  theme_minimal()

plot(p.ratio)


library(gridExtra)

plot_all = grid.arrange(p.rea, p.ne,p.ratio, ncol=1)
plot(plot_all)

ggsave(plot=p.rea,paste("../../Reassortment-Text/Figures/reassortmentRates.pdf", sep=""),width=8, height=4)


