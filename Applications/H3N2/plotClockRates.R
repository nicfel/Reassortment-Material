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

# get all logfiled
log <- list.files(path="./outgroup/", pattern="yh3n2.*norea\\_rep0\\.log", full.names = TRUE)


first = T;
for(i in seq(1, length(log))){
  # set the log file name
  fname1 = log[[i]]
  fname2 = gsub("rep0", "rep1", log[[i]])
  fname3 = gsub("rep0", "rep2", log[[i]])
  
  t.1 = read.table(fname1, header=TRUE, sep="\t")
  t.2 = read.table(fname2, header=TRUE, sep="\t")
  t.3 = read.table(fname3, header=TRUE, sep="\t")
  
  t.1 = t.1[-seq(1,length(t.1$Sample)/5),]
  t.2 = t.2[-seq(1,length(t.2$Sample)/5),]
  t.3 = t.3[-seq(1,length(t.3$Sample)/5),]
  
  t.norea = rbind(t.1,t.2,t.3)
  
  t.1 = read.table(gsub("norea", "", fname1), header=TRUE, sep="\t")
  t.2 = read.table(gsub("norea", "", fname2), header=TRUE, sep="\t")
  t.3 = read.table(gsub("norea", "", fname3), header=TRUE, sep="\t")
  
  t.1 = t.1[-seq(1,length(t.1$Sample)/5),]
  t.2 = t.2[-seq(1,length(t.2$Sample)/5),]
  t.3 = t.3[-seq(1,length(t.3$Sample)/5),]
  
  t = rbind(t.1,t.2,t.3)
  
  hpd.norea = HPDinterval(as.mcmc(t.norea$clockRate.c))
  hpd = HPDinterval(as.mcmc(t$clockRate.c))

    # get the year of the outgroup
  tmp = strsplit(fname1, split="_")[[1]][[2]]
  tmp = strsplit(tmp, split="norea")[[1]][[1]]

  new.clock = data.frame(rate.lower=hpd.norea[1,"lower"], 
                                rate.upper=hpd.norea[1,"upper"],
                                rate.mean=mean(t.norea$clockRate.c),
                         method="independent segments",
                         year=as.numeric(tmp)-0.1)
  
  new.clock = rbind(new.clock, data.frame(rate.lower=hpd[1,"lower"], 
                         rate.upper=hpd[1,"upper"],
                         rate.mean=mean(t$clockRate.c),
                         method="with reassortment", 
                         year=as.numeric(tmp)+0.1)
                    )
  

  if (first){
    clock = new.clock
    first = F;
  }else{
    clock = rbind(clock, new.clock)
  }
  
}



p.events <- ggplot(clock) +
  geom_linerange(aes(x=year, ymin=rate.lower, ymax=rate.upper, color=method), size=2) +
  xlab("year of the outgroup") +
  ylab("substitution rate per site and year") +
  theme_minimal() +
  scale_x_continuous(breaks=seq(2006,2016,2))+
  scale_color_brewer(type='qual')
plot(p.events)
ggsave(plot=p.events,paste("../../../Reassortment-Text/Figures/precision/clockrateDifferences.pdf", sep=""),width=7, height=3)


