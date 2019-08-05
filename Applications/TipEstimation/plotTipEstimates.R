######################################################
######################################################
# Here the inferred mean coalescent and migration
# rate ratios are plotted
######################################################
######################################################
library(ggplot2)
library(coda)
library(grid)
library(gridExtra)
library("colorblindr")

# clear workspace
rm(list = ls())

# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

first = T

# get the names of all output files of the first replicate
log <- list.files(path="./out/", pattern="h3n2.*norea\\_rep.*\\.log", full.names = TRUE)

for (i in seq(1,length(log))){
  print(i)
  fname1 = log[[i]]
  fname2 = gsub("norea", "", log[[i]])
  
  if (file.exists(fname2)){
    # read in the log file
    t.norea <- read.table(fname1, header=TRUE, sep="\t", check.names=FALSE)
  
    # read in the log file
    t <- read.table(fname2, header=TRUE, sep="\t", check.names=FALSE)
    
    # check how the tip taht was estimated is called
    all_labels = labels(t)
    for (j in seq(1,length(all_labels[[2]]))){
      if (startsWith(all_labels[[2]][[j]], "height")){
        heightname = all_labels[[2]][[j]]
      }
    }
    
    # get the true value (the starting value)
    true_value = t[1,heightname]
    if (abs(true_value-t.norea[1,heightname])>0.001){
      print(true_value)
      print(t.norea[1,heightname])
      print("error in the starting value")
    }
  
      # take a 10% burnin
    t <- t[-seq(1,ceiling(length(t$posterior)/10)), ]
    # take a 10% burnin
    t.norea <- t.norea[-seq(1,ceiling(length(t.norea$posterior)/10)), ]
  
  
    # calculate ess values
    ess <- effectiveSize(t)
    ess.norea <- effectiveSize(t.norea)
    
    # non estimated things will have ess 0 
    ess[ess<1]=1000
    ess.norea[ess.norea<1] = 1000
    
    if (min(ess[2:length(ess)])>50 && min(ess.norea[2:length(ess.norea)])>50){
      # get the name of the virus
      virusname = gsub("./out//", "", strsplit(fname2, split="_")[[1]][[1]])
      hpd = HPDinterval(as.mcmc(t[,heightname]))
      hpd.norea = HPDinterval(as.mcmc(t.norea[,heightname]))
      
      cov = hpd[1,'lower']<= true_value & hpd[1,'upper']>= true_value
      cov.norea = hpd.norea[1,'lower']<= true_value & hpd.norea[1,'upper']>= true_value
      # combine with the other replicates
      new.tip = data.frame(true=true_value, 
                          estimated=mean(t[,heightname]), 
                          upper=hpd[1,'upper'],
                          lower=hpd[1,'lower'], 
                          estimated.norea=mean(t.norea[,heightname]), 
                          upper.norea=hpd.norea[1,'upper'],
                          lower.norea=hpd.norea[1,'lower'],
                          run=i,
                          cov =cov,
                          cov.norea=cov.norea,
                          virus=virusname
                          )
      
      hpd = HPDinterval(as.mcmc(t[,"clockRate.c"]))
      hpd.norea = HPDinterval(as.mcmc(t.norea[,"clockRate.c"]))
      
      new.clock = data.frame(estimated=mean(t[,"clockRate.c"]), 
                             upper=hpd[1,'upper'],
                             lower=hpd[1,'lower'], 
                             estimated.norea=mean(t.norea[,"clockRate.c"]), 
                             upper.norea=hpd.norea[1,'upper'],
                             lower.norea=hpd.norea[1,'lower'],
                             run=i,
                             virus=virusname
                             )
      
      hpd = HPDinterval(as.mcmc(t[,"popSize.t"]))
      hpd.norea = HPDinterval(as.mcmc(t.norea[,"popSize.t"]))
      
      new.ne = data.frame(estimated=mean(t[,"popSize.t"]), 
                             upper=hpd[1,'upper'],
                             lower=hpd[1,'lower'], 
                             estimated.norea=mean(t.norea[,"popSize.t"]), 
                             upper.norea=hpd.norea[1,'upper'],
                             lower.norea=hpd.norea[1,'lower'],
                             run=i,
                             virus=virusname
      )
      
  
      if (first){
        tip = new.tip
        clock = new.clock
        ne = new.ne
        first = F
      }else{
        tip = rbind(tip, new.tip)
        clock = rbind(clock, new.clock)
        ne = rbind(ne, new.ne)
      }
    }else{
      print("ess low")
    }
  }
}

# remove tip estimates where the presumably true value and the estimate for both methods differe by more than on year
# tip = tip[-which(abs(tip$estimated-tip$true)>1 & abs(tip$estimated.norea-tip$true)>1),]



p.rea <- ggplot(tip)+
  geom_hline(yintercept = 0, color="red")+
  geom_errorbar(aes(x=run, ymin=lower-true, ymax=upper-true), colour="grey", width=0.01) +
  geom_point(aes(x=run, y=estimated-true, color=virus), size=2) +
  scale_y_continuous(limits=c(-4,4))+
  scale_color_OkabeIto()+
  ylab("difference between estimated and true tip date")+
  xlab("run number")+
  ggtitle(paste("coalescent with reassortment ",
                "(coverage=", sprintf("%.2f", mean(tip$cov)), "",
                ")", sep=""))+
  annotate("text", x=250, y=2, label=paste("coverage for 2 year interval =", sprintf("%.2f", mean(tip[which(tip$virus=="h3n2new"),]$cov)) ))+
  annotate("text", x=750, y=2, label=paste("coverage for 10 year interval =", sprintf("%.2f", mean(tip[which(tip$virus=="h3n2"),]$cov)) ))+
  theme_light() +
  theme(legend.position = "none")
plot(p.rea)

p.norea <- ggplot(tip)+
  geom_hline(yintercept = 0, color="red")+
  geom_errorbar(aes(x=run, ymin=lower.norea-true, ymax=upper.norea-true), colour="grey", width=0.01) +
  geom_point(aes(x=run, y=estimated.norea-true, color=virus), size=2) +
  scale_y_continuous(limits=c(-4,4))+
  scale_color_OkabeIto()+
  ylab("difference between estimated and true tip date")+
  xlab("run number")+
  ggtitle(paste("coalescent with independent segments ",
                "(coverage=", sprintf("%.2f", mean(tip$cov.norea)), "",
                ")", sep=""))+
  annotate("text", x=250, y=2, label=paste("coverage for 2 year interval =", sprintf("%.2f", mean(tip[which(tip$virus=="h3n2new"),]$cov.norea)) ))+
  annotate("text", x=750, y=2, label=paste("coverage for 10 year interval =", sprintf("%.2f", mean(tip[which(tip$virus=="h3n2"),]$cov.norea)) ))+
  theme_light() +
  theme(legend.position = "none")
plot(p.norea)

plot.both <- do.call("grid.arrange", c(list(p.rea, p.norea), ncol=1))
ggsave(plot=plot.both, "../../../Reassortment-Text/Figures/tipest.pdf",width=15, height=7)


p.clock.comp <- ggplot(clock)+
  geom_abline(yintercept = 0, color="red")+
  geom_errorbar(aes(x=estimated, ymin=lower.norea, ymax=upper.norea), colour="grey") +
  geom_errorbarh(aes(y=estimated.norea, xmin=lower, xmax=upper), colour="grey") +
  geom_point(aes(x=estimated, y=estimated.norea, color=virus), size=2) +
  scale_x_continuous(limits=c(0.001,0.004))+
  scale_y_continuous(limits=c(0.001,0.004))+
  xlab("estimated with reassortment")+
  ylab("estimated without reassortment")+
  theme_light()
plot(p.clock.comp)
ggsave(plot=p.clock.comp, "../../../Reassortment-Text/Figures/precision/clockcomp.pdf",width=15, height=7)


p.Ne.comp <- ggplot(ne)+
  geom_abline(yintercept = 0, color="red")+
  geom_errorbar(aes(x=estimated, ymin=lower.norea, ymax=upper.norea), colour="grey") +
  geom_errorbarh(aes(y=estimated.norea, xmin=lower, xmax=upper), colour="grey") +
  geom_point(aes(x=estimated, y=estimated.norea, color=virus), size=2) +
  xlab("estimated with reassortment")+
  ylab("estimated without reassortment")+
  scale_x_continuous(limits=c(0,12))+
  scale_y_continuous(limits=c(0,12))+
  theme_light()
plot(p.Ne.comp)
ggsave(plot=p.Ne.comp, "../../../Reassortment-Text/Figures/precision/Necomp.pdf",width=15, height=7)


