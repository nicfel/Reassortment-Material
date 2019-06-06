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
log <- list.files(path="./out/", pattern="h3n2newnorea\\_rep.*\\.log", full.names = TRUE)

for (i in seq(1,length(log))){
  print(i)
  fname1 = log[[i]]
  fname2 = gsub("norea", "", log[[i]])
  
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

    if (first){
      tip = new.tip
      first = F
    }else{
      tip = rbind(tip, new.tip)
    }
  }else{
    print("ess low")
  }
}

# remove tip estimates where the presumably true value and the estimate for both methods differe by more than on year
tip = tip[-which(abs(tip$estimated-tip$true)>1 & abs(tip$estimated.norea-tip$true)>1),]



p.rea <- ggplot(tip)+
  geom_hline(yintercept = 0, color="red")+
  geom_errorbar(aes(x=run, ymin=lower-true, ymax=upper-true), colour="grey", width=0.01) +
  geom_point(aes(x=run, y=estimated-true, color=virus), size=2) +
  scale_y_continuous(limits=c(-4,4))+
  scale_color_OkabeIto()+
  ylab("difference between estimated and true tip date")+
  xlab("run number")+
  ggtitle(paste("coalescent with reassortment ",
                "(coverage: overall=", sprintf("%.2f", mean(tip$cov)), ", ",
                "p09H1N1=", sprintf("%.2f", mean(tip[which(tip$virus=="h1n1pdm"),]$cov)), ", ",
                "seasH1N1=", sprintf("%.2f", mean(tip[which(tip$virus=="h1n1sea"),]$cov)), ", ",
                "H3N2=", sprintf("%.2f", mean(tip[which(tip$virus=="h3n2"),]$cov)), ", ",
                "InfB=", sprintf("%.2f", mean(tip[which(tip$virus=="infB"),]$cov)), 
                ")", sep=""))+
  theme_light()
plot(p.rea)

p.norea <- ggplot(tip)+
  geom_hline(yintercept = 0, color="red")+
  geom_errorbar(aes(x=run, ymin=lower.norea-true, ymax=upper.norea-true), colour="grey", width=0.01) +
  geom_point(aes(x=run, y=estimated.norea-true, color=virus), size=2) +
  scale_y_continuous(limits=c(-4,4))+
  scale_color_OkabeIto()+
  ylab("difference between estimated and true tip date")+
  xlab("run number")+
  ggtitle(paste("same coalescent process, independent trees ",                
                "(coverage: overall=", sprintf("%.2f", mean(tip$cov.norea)), ", ",
                "p09H1N1=", sprintf("%.2f", mean(tip[which(tip$virus=="h1n1pdm"),]$cov.norea)), ", ",
                "seasH1N1=", sprintf("%.2f", mean(tip[which(tip$virus=="h1n1sea"),]$cov.norea)), ", ",
                "H3N2=", sprintf("%.2f", mean(tip[which(tip$virus=="h3n2"),]$cov.norea)), ", ",
                "InfB=", sprintf("%.2f", mean(tip[which(tip$virus=="infB"),]$cov.norea)), 
                ")", sep=""))+
  
  theme_light()
plot(p.norea)

plot.both <- do.call("grid.arrange", c(list(p.rea, p.norea), ncol=1))
ggsave(plot=plot.both, "../../../Reassortment-Text/Figures/precision/tipest.pdf",width=15, height=7)

  