######################################################
######################################################
# combine the gene tree runs and run the mcc trees
######################################################
######################################################
# clear workspace
library(ggplot2)
library(grid)
library(gridExtra)
library(colorblindr)
rm(list = ls())

# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

heights <- list.files(path=paste("./trees/heights",sep=""), pattern="*\\.txt$", full.names = TRUE)
t = read.table("./trees/Heights.csv", header=T, sep=",", na.strings="")

t$mid1 = (t$lower1+t$upper1)/2
t$mid2 = (t$lower2+t$upper2)/2

t$ratio = (t$upper2-t$lower2)/(t$upper1-t$lower1)

t = t[-which(t$post1<0.25 | t$post2<0.25),]

p <- ggplot(data=t) +
  geom_histogram(aes(ratio),color="black", fill="white") +
  geom_vline(xintercept=median(t$ratio), color="red")+
  scale_x_log10(limits=c(0.1,10)) +
  annotate("text", x=1.5,y=200,label=sprintf("median=%.2f",median(t$ratio)))+
  theme_minimal()+
  xlab("reduction in uncertainty") +
  theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
  )


# read in the log files as well to compare clock rates and effective population sizes
t.rea = read.table("./combined/h3n2.log", header=T, sep="\t", na.strings="")
t.norea = read.table("./combined/h3n2norea.log", header=T, sep="\t", na.strings="")

p.ne <- ggplot() +
  geom_density(data=t.rea,aes(popSize.t, fill=" with reassortment")) +
  geom_density(data=t.norea,aes(popSize.t, fill=" independent segments")) +
  theme_minimal()+
  xlab("effective population size")+
  scale_fill_brewer() +
  theme(legend.position="none",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
  )

plot(p.ne)

p.clock <- ggplot() +
  geom_density(data=t.rea, aes(clockRate.c, fill=" with reassortment")) +
  geom_density(data=t.norea, aes(clockRate.c, fill=" independent segments")) +
  theme_minimal()+
  xlab("clock rate")+
  scale_fill_brewer() +
  theme(legend.title = element_blank(),
        legend.position="top",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
        )

plot(p.clock)

                             
  
g_legend <- function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  return(legend)} 

legend.plot <- p.clock + theme(legend.position="right")
legend <- g_legend(legend.plot) 

p.clock = p.clock +theme(legend.position="none")


plot.all = grid.arrange(p, p.ne, legend, p.clock, ncol=2)


ggsave(plot=plot.all,paste("../../../Reassortment-Text/Figures/precision_bias.pdf", sep=""),width=6, height=6)


# plot.mean = do.call("grid.arrange",c(p.mean, ncol=8))
# plot(p.mean)

