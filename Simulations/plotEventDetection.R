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

# read in true 
events = read.table("event_detection.csv", header=TRUE, sep=",")

p.true <- ggplot(events,aes(x=max, y=max-min,z=detection_probability))+
  stat_contour(geom = "polygon", aes(fill = ..level..)) +
  geom_tile(aes(fill = detection_probability)) +
  stat_contour(bins = 15) +
  theme_minimal() 
  # scale_x_log10() 
plot(p.true)
p.true <- ggplot(events,aes(x=max, y=max-min,z=detection_probability))+
  geom_tile(aes(fill = detection_probability)) +
  theme_minimal() 
# scale_x_log10() 
plot(p.true)


# read in false events
noevent = read.table("noevent_detection.csv", header=TRUE, sep=",")

p.false <- ggplot(noevent)+
  geom_point(aes(x=max, y=difference,color=detection_probability), size=1) + 
  theme_minimal()
# plot(p.false)


# 
# 
# # get the names of all output files of the first replicate
# log <- list.files(path="./out", pattern="*rep1.log", full.names = TRUE)
# 
# # use the matlab standard colors to plot
# col0 <- rgb(red=0.0, green=0.4470,blue=0.7410)
# col1 <- rgb(red=0.8500, green=0.3250,blue=0.0980)
# col2 <- rgb(red=0.9290, green=0.6940,blue=0.1250)
# col4 <- rgb(red=0.4660, green=0.6740,blue=0.1880)
# col3 <- rgb(red=0.3010, green=0.7450,blue=0.9330)
# 
# first = T
# c = 0;
# for (i in seq(1,length(log))){
#   # read in the log file
#   t.1 <- read.table(log[[i]], header=TRUE, sep="\t")
#   t.2 <- read.table(gsub("rep1","rep2",log[[i]]), header=TRUE, sep="\t")
#   t.3 <- read.table(gsub("rep1","rep3",log[[i]]), header=TRUE, sep="\t")
#   
#   # take a 10% burnin
#   t.1 <- t.1[-seq(1,ceiling(length(t.1$posterior)/5)), ]
#   t.2 <- t.3[-seq(1,ceiling(length(t.2$posterior)/5)), ]
#   t.3 <- t.3[-seq(1,ceiling(length(t.3$posterior)/5)), ]
#   
#   # combine the replictes
#   t = rbind(t.1,t.2,t.3)
#   # calculate ess values
#   ess <- effectiveSize(t)
#   
#   post_ess = as.numeric(ess["reassortmentRate"])
#   
#   if (post_ess>200){
#     
#     # find the correct run number
#     tmp = strsplit(log[[i]], '_')[[1]][[2]]
#     used.rates = true.rates[which(true.rates$run==as.numeric(tmp)),];
#   
#     # combine with the other replicates
#     new.rea = data.frame(true=used.rates$reassortment, estimated=mean(t$reassortmentRate), 
#                     upper=quantile(t$reassortmentRate,0.975), lower=quantile(t$reassortmentRate,0.025) )
#     new.Ne = data.frame(true=used.rates$Ne, estimated=mean(t$popSize), 
#                     upper=quantile(t$popSize,0.975), lower=quantile(t$popSize,0.025) )
#     
#     
#     if (first){
#       rea = new.rea
#       Ne = new.Ne
#       first = F
#     }else{
#       rea = rbind(rea, new.rea)
#       Ne = rbind(Ne, new.Ne)
#     }
#   }else{
#     c = c+1
#     print(c)
#   }
# }
# 
# p.Ne <- ggplot(Ne)+
#   geom_abline(intercept = 0, color="red")+
#   geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
#   geom_point(aes(x=true, y=estimated), size=2) + 
#   theme_minimal() +
#   scale_y_log10() +
#   scale_x_log10() 
# p.rea <- ggplot(rea)+
#   geom_abline(intercept = 0, color="red", linetype="dashed")+
#   geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
#   geom_point(aes(x=true, y=estimated), size=2) + 
#   theme_minimal() +
#   scale_y_log10() +
#   scale_x_log10() 
# 
#   # scale_y_log10(limits=c(0.01,1.0)) +
#   # scale_x_log10(limits=c(0.01,1.0)) 
# 
# plot(p.Ne)
# plot(p.rea)
# 
# ggsave(plot=p.Ne,paste("../../Reassortment-Text/Figures/Ne_sim.pdf", sep=""),width=6, height=5)
# ggsave(plot=p.rea,paste("../../Reassortment-Text/Figures/rho_sim.pdf", sep=""),width=6, height=5)
