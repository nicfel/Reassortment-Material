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


gen_information = c("low", "high")

for (g in seq(1,length(gen_information))){
  # get the names of all output files of the first replicate
  log <- list.files(path="./out/", pattern=paste("inf_", gen_information[[g]], ".*rep1.log", sep=""), full.names = TRUE)
  
  first = T
  c = 0;
  for (i in seq(1,length(log))){
    # read in the log file
    t.1 <- read.table(log[[i]], header=TRUE, sep="\t")
    t.2 <- read.table(gsub("rep1","rep2",log[[i]]), header=TRUE, sep="\t")
    t.3 <- read.table(gsub("rep1","rep3",log[[i]]), header=TRUE, sep="\t")
    
    # take a 10% burnin
    t.1 <- t.1[-seq(1,ceiling(length(t.1$posterior)/5)), ]
    t.2 <- t.3[-seq(1,ceiling(length(t.2$posterior)/5)), ]
    t.3 <- t.3[-seq(1,ceiling(length(t.3$posterior)/5)), ]
    
    # combine the replictes
    t = rbind(t.1,t.2,t.3)
    # calculate ess values
    ess <- effectiveSize(t)
    
    post_ess = as.numeric(ess["posterior"])
    
    if (post_ess>200){
      
      # find the correct run number
      tmp = strsplit(log[[i]], '_')[[1]][[3]]
      used.rates = true.rates[which(true.rates$run==as.numeric(tmp)),];
    
      # combine with the other replicates
      new.rea = data.frame(true=used.rates$reassortment, estimated=mean(t$reassortmentRate), 
                      upper=quantile(t$reassortmentRate,0.975), lower=quantile(t$reassortmentRate,0.025) )
      new.Ne = data.frame(true=used.rates$Ne, estimated=mean(t$popSize), 
                      upper=quantile(t$popSize,0.975), lower=quantile(t$popSize,0.025))
      
      
      if (first){
        rea = new.rea
        Ne = new.Ne
        first = F
      }else{
        rea = rbind(rea, new.rea)
        Ne = rbind(Ne, new.Ne)
      }
    }else{
      c = c+1
      print(c)
    }
  }
  
  p.Ne <- ggplot(Ne)+
    geom_abline(intercept = 0, color="red")+
    geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
    geom_point(aes(x=true, y=estimated), size=2) + 
    theme_minimal() +
    scale_y_log10() +
    scale_x_log10() 
  p.rea <- ggplot(rea)+
    geom_abline(intercept = 0, color="red", linetype="dashed")+
    geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
    geom_point(aes(x=true, y=estimated), size=2) + 
    theme_minimal() +
    scale_y_log10() +
    scale_x_log10() 
  
    # scale_y_log10(limits=c(0.01,1.0)) +
    # scale_x_log10(limits=c(0.01,1.0)) 
  
  plot(p.Ne)
  plot(p.rea)
  ggsave(plot=p.Ne,paste("../../Reassortment-Text/Figures/Ne_sim_", gen_information[[g]], ".pdf", sep=""),width=6, height=5)
  ggsave(plot=p.rea,paste("../../Reassortment-Text/Figures/rho_sim_", gen_information[[g]], ".pdf", sep=""),width=6, height=5)
  
}

