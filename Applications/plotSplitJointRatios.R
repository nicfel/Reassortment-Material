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

viruses = c("H1N1", "H1N1wop", "H3N2", "InfB")
nr_segments = c(8,8,8,8)

first = T;
for(i in seq(1, length(viruses))){
  # get the names run log files
  log <- list.files(path=paste("./",viruses[[i]] ,"/out", sep=""), pattern="*Stats.log", full.names = TRUE)
  # read them in
  for (j in seq(1, length(log))){
    t <- read.table(log[j], header=TRUE, sep="\t")
    # take a 10 % burn in
    t <- t[-seq(1,ceiling(length(t$Sample)/2)), ]
    
    first_reassort = T
    
    for (a in seq(1,nr_segments[[i]]-1)){
      for (b in seq(a+1,nr_segments[[i]])){
        jointname = paste("network.jointReassortment.", a-1, "_", b-1, sep="")
        splitname = paste("network.splitReassortment.", a-1, "_", b-1, sep="")
        new.new.reassortment = data.frame(segment0=a-1,segment1=b-1, ratio=median(log10(t[,jointname]/t[,splitname])))
        new.new.reassortment = rbind(new.new.reassortment, data.frame(segment0=b-1,segment1=a-1, ratio=median(log10(t[,jointname]/t[,splitname]))))
        
        if(a==1&&b==2){
          new.reassortment =new.new.reassortment
        }else{
          new.reassortment = rbind(new.reassortment,new.new.reassortment)
        }
      }
    }
    # normalize the results
    new.reassortment$ratio = new.reassortment$ratio - mean(new.reassortment$ratio)
    
    
    new.reassortment$virus = viruses[[i]]
    
    if (first){
      reassortment = new.reassortment
      first = F;
    }else{
      reassortment = rbind(reassortment, new.reassortment)
    }
  }
}

segments = c("HA", "MP", "NA", "NP", "NS1", "PA", "PB1", "PB2")

p <- ggplot(reassortment) +
  geom_tile(aes(segment0,segment1,fill=ratio))+
  scale_fill_gradientn(colors =  c("#88419d","#b3cde3", rgb(1,1,1), "#fdcc8a","#d7301f"))+
  facet_grid(.~virus)+
  theme_minimal() +
  scale_x_continuous(breaks=seq(0,7), labels=segments)+
  scale_y_continuous(breaks=seq(0,7), labels=segments)
  
plot(p)
  

# p.rea <- ggplot(reassortment) +
#   geom_violin(aes(x=virus, y=rate)) +
#   scale_y_log10() +
#   xlab("") +
#   ylab("reassortment rate") +
#   theme_minimal()
# 
# p.ne <- ggplot(reassortment) +
#   geom_violin(aes(x=virus, y=Neff)) +
#   scale_y_log10() +
#   xlab("") +
#   ylab("effective population size") +
#   theme_minimal()
# 
# p.ratio <- ggplot(reassortment) +
#   geom_violin(aes(x=virus, y=rate/Neff)) +
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
# plot_all = grid.arrange(p.rea, p.ne,p.ratio, ncol=1)
# plot(plot_all)
# 
# ggsave(plot=plot_all,paste("../../Reassortment-Text/Figures/reassortmentRates.pdf", sep=""),width=8, height=8)


