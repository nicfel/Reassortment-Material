######################################################
######################################################
# Here the inferred mean coalescent and migration
# rate ratios are plotted
######################################################
######################################################
library(ggplot2)
library("coda")
library("colorblindr")

# clear workspace
rm(list = ls())

# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)


# viruses = c("H1N1pandemic", "H1N1seasonal", "H3N2", "InfB", "InfC", "InfD")
# 
# virus_labels = c("p2009 like H1N1", "seasonal H1N1", "H3N2", "H5N1", "Influenza B yam + vic", "Influenza C", "Influenza D")
# 
# 
# 
# viruses = c("H1N1pandemic", "H1N1seasonal", "H3N2", "H3N2", "H3N2", "H3N2", "H3N2", "H2N2", "InfB")
# 
# virus_labels = c("p09 like H1N1", "pre 09 H1N1", "H3N2 1980-2020", "H3N2 1980-1990", "H3N2 1990-2000", "H3N2 2000-2010", "H3N2 2010-2020", "H2N2", "Influenza B")
# 
# virusname = c("h1n1pdm", "h1n1sea", "h3n2", "h3n2ancient", "h3n2old","h3n2recent","h3n2new", "h2n2", "infB")


viruses = c("H1N1seasonal","H1N1seasonal","H1N1seasonal",
            "H1N1pandemic","H1N1pandemic","H1N1pandemic",
            "H3N2", "H3N2", "H3N2",
            "H2N2",
            "InfB", "InfB", "InfB")

virus_labels = c("pre 09 H1N1", "pre 09 H1N1",  "pre 09 H1N1", 
                 "p09 like H1N1", "p09 like H1N1",  "p09 like H1N1",  
                 "H3N2", "H3N2", "H3N2",
                 "H2N2",  
                 "Influenza B", "Influenza B", "Influenza B")

virusname = c("h1n1sea_sub1","h1n1sea_sub2","h1n1sea_sub3",
              "h1n1pdm_sub1","h1n1pdm_sub2","h1n1pdm_sub3", 
              "h3n2_sub1", "h3n2_sub2", "h3n2_sub3", 
              "h2n2", 
              "infB_sub1", "infB_sub2", "infB_sub3")

replicate = c("1","2","3",
              "1","2","3", 
              "1", "2", "3", 
              "2", 
              "1", "2", "3")
offset=0.2

xval = c(1-offset,1,1+offset,
         2-offset,2,2+offset,
         3-offset,3,3+offset,
         4,
         5-offset,5,5+offset)




first = T;
for(i in seq(1, length(viruses))){
  # set the log file name
  filename = path=paste("./",viruses[[i]] ,"/combined/", virusname[[i]],".log", sep="")
  
  # read the file in
  t <- read.table(filename, header=TRUE, sep="\t")
  
  # get the HPD's
  hpd.rea = HPDinterval(as.mcmc(t$reassortmentRate))
  hpd.ne = HPDinterval(as.mcmc(t$popSize.t))
  hpd.clock = HPDinterval(as.mcmc(t$clockRate.c))
  hpd.eventsyear = HPDinterval(as.mcmc(t$network.obsReassortmentNodeCount/t$network.obsHeight))
  
  new.reassortment = data.frame(rate.lower=hpd.rea[1,"lower"], rate.upper=hpd.rea[1,"upper"],
                                ne.lower=hpd.ne[1,"lower"], ne.upper=hpd.ne[1,"upper"],
                                clock.lower=hpd.clock[1,"lower"], clock.upper=hpd.clock[1,"upper"],
                                events.lower=hpd.eventsyear[1,"lower"], events.upper=hpd.eventsyear[1,"upper"])
  
  new.reassortment$virus = virus_labels[[i]]
  new.reassortment$replicate = replicate[[i]]
  new.reassortment$xval = xval[[i]]
  
  if (first){
    reassortment = new.reassortment
    first = F;
  }else{
    reassortment = rbind(reassortment, new.reassortment)
  }
  
}



# reassortment$virus = factor(reassortment$virus, levels = virus_labels)


tick_labels = c("p1918 like H1N1",
                "p2009 like H1N1",
                "H3N2",
                "H2N2",  
                "Influenza B")

p.events <- ggplot(reassortment) +
  geom_linerange(aes(x=xval, ymin=events.lower, ymax=events.upper, color=replicate), size=4) +
  xlab("") +
  ylab("number of events per year") +
  theme_minimal() +
  scale_x_continuous(breaks=seq(1,5), label=tick_labels) +
  scale_y_continuous(limits=c(0,8))+
  scale_color_viridis_d()
plot(p.events)
ggsave(plot=p.events,paste("../../Reassortment-Text/Figures/networks/eventsyear.pdf", sep=""),width=7, height=3)


p.rea <- ggplot(reassortment) +
  geom_linerange(aes(x=xval, ymin=rate.lower, ymax=rate.upper, color=replicate), size=4) +
  xlab("") +
  ylab("reassortment rate") +
  theme_minimal() +
  scale_x_continuous(breaks=seq(1,5), label=tick_labels) +
  scale_color_viridis_d()
plot(p.rea)
ggsave(plot=p.rea,paste("../../Reassortment-Text/Figures/networks/reassortmentRates.pdf", sep=""),width=7, height=3)


# p.ne <- ggplot(reassortment) +
#   geom_linerange(aes(x=virus, ymin=ne.lower, ymax=ne.upper), size=10) +
#   scale_y_log10() +
#   xlab("") +
#   ylab("effective population size") +
#   theme_minimal()
# p.clock <- ggplot(reassortment) +
#   geom_violin(aes(x=virus, y=clock)) +
#   scale_y_log10() +
#   xlab("") +
#   ylab("clock rate") +
#   theme_minimal()
# 
# p.ratio <- ggplot(reassortment) +
#   geom_violin(aes(x=virus, y=rate/clock)) +
#   scale_y_log10() +
#   xlab("") +
#   ylab("reassortment rate over Ne") +
#   theme_minimal()
# 
# plot(p.ratio)
# 
# 
# library(gridExtra)
# 
# plot_all = grid.arrange(p.rea, p.clock, p.ne,p.ratio, ncol=2)
# plot(plot_all)
# 
# 
# 
