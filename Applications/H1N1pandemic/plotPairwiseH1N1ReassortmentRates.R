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

segments = c("HA", "MP", "NA", "NP", "NS1", "PA", "PB1", "PB2", "prior")

prior_val = c(0.01)

first = T;
for (a in seq(1,length(segments)-2)){
  for (b in seq(a+1,length(segments)-1)){
    # get the names run log files
    print("i")
    log <- list.files(path="./pairout/", pattern=paste("h1n1pdm_",segments[[a]], "_", segments[[b]], "_rep", "*[012].log",sep=""), full.names = TRUE)
    log <- list.files(path="./pairout/", pattern=paste("h1n1pdm_",segments[[a]], "_", segments[[b]], "_rep", "*[012].log",sep=""), full.names = TRUE)
    print("o")
    for (j in seq(1, length(log))){
      t <- read.table(log[j], header=TRUE, sep="\t")
      # take a 10 % burn in
      t <- t[-seq(1,ceiling(length(t$Sample)/10)), ]
      new.rate = data.frame(rate=t$reassortmentRate, obscount=t$network.obsReassortmentNodeCount,segment1=segments[[a]],segment2=segments[[b]])
      new.rate = rbind(new.rate, data.frame(rate=t$reassortmentRate, obscount=t$network.obsReassortmentNodeCount, segment1=segments[[b]],segment2=segments[[a]]))
      if (first){
        rate=new.rate
        first=F
      }else{
        rate=rbind(rate,new.rate)
      }
    }
  }
}
# # add prior
# for (a in seq(1,length(segments))){
#   for (b in seq(1,1000)){
#     new.rate = data.frame(rate=rlnorm(1,0,4),segment1=segments[[a]],segment2=segments[[9]])
#     rate=rbind(rate,new.rate)
# 
#   }
# }
# 

# set the levels 
rate$segment2 <- factor(rate$segment2, levels = c("HA", "NA", "MP", "NP", "NS1", "PA", "PB1", "PB2"))
rate$segment1 <- factor(rate$segment1, levels = c("HA", "NA", "MP", "NP", "NS1", "PA", "PB1", "PB2"))


p.rea <- ggplot(rate) +
  geom_violin(aes(x=segment2, y=rate,fill=segment2),colour="black") +
  scale_y_log10(limits=c(0.01,1)) +
  facet_wrap(~segment1)+
  xlab("") +
  ylab("reassortment rate") +
  scale_fill_manual(name="",values=c('#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02','#a6761d','#666666', '#ffffff')) +
  theme(legend.position="top")
require(lemon)
p.rea = reposition_legend(p.rea, 'center', panel='panel-3-3')

plot(p.rea)
ggsave(plot=p.rea,paste("../../../Reassortment-Text/Figures/pairwise/h1n1pdm_pairs", ".pdf" , sep=""),width=10, height=6)


# rate$segment2 <- factor(rate$segment2, levels = c("HA", "NA", "MP", "NP", "NS1", "PA", "PB1", "PB2"))
# 
# p.rea.count <- ggplot(rate) +
#   geom_violin(aes(x=segment2, y=obscount,fill=segment2),colour="black") +
#   scale_y_continuous(limits=c(0,20)) +
#   facet_wrap(~segment1)+
#   xlab("") +
#   ylab("reassortment rate") +
#   scale_fill_manual(name="",values=c('#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02','#a6761d','#666666')) +
#   theme(legend.position="top")
# require(lemon)
# p.rea.count = reposition_legend(p.rea.count, 'center', panel='panel-3-3')
# 
# plot(p.rea.count)
# ggsave(plot=p.rea.count,paste("../../../Reassortment-Text/Figures/pairwise/h1n1_count", ".pdf" , sep=""),width=10, height=6)



# p.ind <- ggplot(rate) +
#   geom_violin(aes(x=segment2, y=rate,fill=segment2),colour="black") +
#   scale_y_log10(limits=c(0.1,1)) +
#   xlab("") +
#   ylab("reassortment rate") +
#   scale_fill_manual(name="",values=c('#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02','#a6761d','#666666')) +
#   theme(legend.position="none")
# 
# plot(p.ind)
# ggsave(plot=p.ind,paste("../../../Reassortment-Text/Figures/pairwise/h3n2_ind", ".pdf" , sep=""),width=5, height=3)
#   
# 
# 
# 
# library(gridExtra)
# 
# plot_all = grid.arrange(p.rea, p.ne,p.ratio, ncol=1)
# plot(plot_all)
# 
# ggsave(plot=p.rea,paste("../../Reassortment-Text/Figures/reassortmentRates.pdf", sep=""),width=8, height=4)


