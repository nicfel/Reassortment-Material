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
library(coda)

# clear workspace
rm(list = ls())

# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

logs <- list.files(path=paste("./shortout/combined/",sep=""), pattern=".*\\.log$", full.names = TRUE)

first = T;
for(i in seq(1, length(logs))){
  # read in the log file
  t <- read.table(logs[[i]], header=T, sep="\t")
  # get the year
  tmp = strsplit(logs[[i]], split="_")[[1]][[2]]
  tmp = strsplit(tmp, split="\\.")[[1]][[1]]
  # compute HPD
  hpd = HPDinterval(as.mcmc(t$reassortmentRate))
  new.rate = data.frame(year=as.numeric(tmp), lower=hpd[1,"lower"], upper=hpd[1,"upper"], mean=mean(t$reassortmentRate))

  if (first){
    rate = new.rate
    first = F;
  }else{
    rate = rbind(rate, new.rate)
  }
}

trunktime = c("2","4","6")
subset = c("1","2","3")

# read in the trunk reassortment rates
lsPlots = list()
lower.vals.fit = list()
upper.vals.fit = list()

lower.vals.unfit = list()
upper.vals.unfit = list()

c=1
grid.a.b <- expand.grid(seq(1,length(subset)), seq(1,length(trunktime)))

lsPlots <- lapply(seq_len(nrow(grid.a.b)), function(i) {
      a <- grid.a.b[i, "Var1"]
      b <- grid.a.b[i, "Var2"]
    t <- read.table(paste("./combined/h3n2_sub",subset[[a]], ".", trunktime[[b]], "trunk.txt", sep=""), header=F, sep="\t")
    
    fit.rea.rate = t$V1/t$V3
    unfit.rea.rate = t$V2/t$V4
    
    hpd.fit = HPDinterval(as.mcmc(fit.rea.rate))
    hpd.unfit = HPDinterval(as.mcmc(unfit.rea.rate))
    
    lower.vals.fit[[c]] = hpd.fit[1,"lower"]
    upper.vals.fit[[c]] = hpd.fit[1,"upper"]
    
    lower.vals.unfit[[c]] = hpd.unfit[1,"lower"]
    upper.vals.unfit[[c]] = hpd.unfit[1,"upper"]
    
    
    p.sub = ggplot(rate) +
      geom_ribbon(aes(x=year, 
                      ymin=rep(lower.vals.fit[[eval(c)]] ,length(rate$year)), 
                      ymax=rep(upper.vals.fit[[eval(c)]] ,length(rate$year)), fill="fit") )+
      
      geom_ribbon(aes(x=year, 
                      ymin=rep(lower.vals.unfit[[c]],length(rate$year)), 
                      ymax=rep(upper.vals.unfit[[c]],length(rate$year)), fill="unfit") )+
      
      scale_fill_OkabeIto()+
      
      
      geom_linerange(aes(x=year,ymin=lower, ymax=upper))+
      geom_point(aes(x=year,y=mean)) +
      scale_x_continuous(limits=c(2000,2016)) + # something odd is happening with 2018
      theme_minimal() +
      theme(legend.position = "none") +
      ylab("reassortment rate per lineage and year") +
      ggtitle(paste("subset", subset[[a]], "tip distance of fit edge = ", trunktime[[b]], "years"))
    p.sub
})

plot.all = do.call("grid.arrange",c(lsPlots, ncol=3))


ggsave(plot=plot.all, paste("../../../Reassortment-Text/Figures/short_term_reassortment.pdf", sep=""),width=15, height=10)
