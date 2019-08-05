######################################################
######################################################
# Here the inferred mean coalescent and migration
# rate ratios are plotted
######################################################
######################################################
library(ggplot2)
library(coda)


# clear workspace
rm(list = ls())

# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

# read in the true rates
true.rates <- read.table("rates.csv", header=TRUE, sep=",")
first = T

# get the names of all output files of the first replicate
log <- list.files(path="./out/", pattern=paste("inf_low", ".*rep1.log", sep=""), full.names = TRUE)

for (i in seq(1,length(log))){
  print(i)
  fname1 = log[[i]]
  fname2 = gsub("low", "high", log[[i]])
  
  # read in the log file
  t.1 <- read.table(fname1, header=TRUE, sep="\t")
  t.2 <- read.table(gsub("rep1","rep2",fname1), header=TRUE, sep="\t")
  t.3 <- read.table(gsub("rep1","rep3",fname1), header=TRUE, sep="\t")
  
  # read in the log file
  t.1.high <- read.table(fname2, header=TRUE, sep="\t")
  t.2.high <- read.table(gsub("rep1","rep2",fname2), header=TRUE, sep="\t")
  t.3.high <- read.table(gsub("rep1","rep3",fname2), header=TRUE, sep="\t")
  
  # take a 10% burnin
  t.1 <- t.1[-seq(1,ceiling(length(t.1$posterior)/5)), ]
  t.2 <- t.3[-seq(1,ceiling(length(t.2$posterior)/5)), ]
  t.3 <- t.3[-seq(1,ceiling(length(t.3$posterior)/5)), ]
  # take a 10% burnin
  t.1.high <- t.1.high[-seq(1,ceiling(length(t.1.high$posterior)/5)), ]
  t.2.high <- t.2.high[-seq(1,ceiling(length(t.2.high$posterior)/5)), ]
  t.3.high <- t.3.high[-seq(1,ceiling(length(t.3.high$posterior)/5)), ]
  
  
  # combine the replictes
  t = rbind(t.1,t.2,t.3)
  t.high = rbind(t.1.high, t.2.high, t.3.high)
  
  # calculate ess values
  ess <- effectiveSize(t)
  ess.high <- effectiveSize(t.high)
  
  if (min(ess[2:length(ess)])>50 && min(ess[2:length(ess.high)])>50){
    
    # find the correct run number
    tmp = strsplit(log[[i]], '_')[[1]][[3]]
    used.rates = true.rates[which(true.rates$run==as.numeric(tmp)),];
  
    # combine with the other replicates
    new.rea = data.frame(true=used.rates$reassortment, 
                        estimated=mean(t$reassortmentRate), 
                        upper=quantile(t$reassortmentRate,0.975),
                        lower=quantile(t$reassortmentRate,0.025), 
                        estimated.high=mean(t.high$reassortmentRate), 
                        upper.high=quantile(t.high$reassortmentRate,0.975),
                        lower.high=quantile(t.high$reassortmentRate,0.025)
                        )
    new.Ne = data.frame(true=used.rates$Ne, 
                        estimated=mean(t$popSize),
                        upper=quantile(t$popSize,0.975), 
                        lower=quantile(t$popSize,0.025), 
                        estimated.high=mean(t.high$popSize),
                        upper.high=quantile(t.high$popSize,0.975), 
                        lower.high=quantile(t.high$popSize,0.025)
                        )
    
    
    if (first){
      rea = new.rea
      Ne = new.Ne
      first = F
    }else{
      rea = rbind(rea, new.rea)
      Ne = rbind(Ne, new.Ne)
    }
  }else{
    print("ess low")
  }
}


p.Ne.comp <- ggplot(Ne)+
  geom_abline(intercept = 0, color="red")+
  geom_errorbar(aes(x=estimated, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_errorbarh(aes(y=estimated.high, xmin=lower.high, xmax=upper.high), colour="grey", width=0.01) +
  geom_point(aes(x=estimated, y=estimated.high), size=2) + 
  theme_light() +
  scale_y_log10() +
  scale_x_log10() +
  xlab("estimated using low evolutionary rate")+
  ylab("estimated using high evolutionary rate")+
  ggtitle("effective population sizes")
plot(p.Ne.comp)

p.rea.comp <- ggplot(rea)+
  geom_abline(intercept = 0, color="red")+
  geom_errorbar(aes(x=estimated, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_errorbarh(aes(y=estimated.high, xmin=lower.high, xmax=upper.high), colour="grey", width=0.01) +
  geom_point(aes(x=estimated, y=estimated.high), size=2) + 
  theme_light() +
  scale_y_log10() +
  scale_x_log10() +
  xlab("estimated using low evolutionary rate")+
  ylab("estimated using high evolutionary rate")+
  ggtitle("reassortment rates")
plot(p.rea.comp)

ggsave(plot=p.Ne.comp,paste("../../../Reassortment-Text/Figures/simulation/Ne_comp.pdf", sep=""),width=4, height=4)
ggsave(plot=p.rea.comp,paste("../../../Reassortment-Text/Figures/simulation/rho_comp.pdf", sep=""),width=4, height=4)


p.Ne <- ggplot(Ne)+
  geom_abline(intercept = 0, color="red")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated), size=2) + 
  theme_minimal() +
  scale_y_log10() +
  scale_x_log10() +
  xlab("true effective population size")+
  ylab("estimated effective population size")
  

p.rea <- ggplot(rea)+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated), size=2) + 
  theme_minimal() +
  scale_y_log10() +
  scale_x_log10() +
  xlab("true reassortment rate")+
  ylab("estimated reassortment rate")


# scale_y_log10(limits=c(0.01,1.0)) +
# scale_x_log10(limits=c(0.01,1.0)) 

plot(p.Ne)
plot(p.rea)
ggsave(plot=p.Ne,paste("../../../Reassortment-Text/Figures/simulation/Ne_sim.pdf", sep=""),width=4, height=4)
ggsave(plot=p.rea,paste("../../../Reassortment-Text/Figures/simulation/rho_sim.pdf", sep=""),width=4, height=4)



