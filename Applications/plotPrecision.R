######################################################
######################################################
# Here the inferred mean coalescent and migration
# rate ratios are plotted
######################################################
######################################################
library(ggplot2)
library("coda")
library("colorblindr")
library(grid)
library(gridExtra)

# clear workspace
rm(list = ls())

# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

viruses = c("H1N1seasonal","H1N1seasonal","H1N1seasonal",
            "H1N1pandemic","H1N1pandemic","H1N1pandemic",
            "H3N2", "H3N2", "H3N2",
            "InfB", "InfB", "InfB")

virus_labels = c("p1918 like H1N1", "p1918 like H1N1",  "p1918 like H1N1", 
                 "p2009 like H1N1", "p2009 like H1N1",  "p2009 like H1N1",  
                 "H3N2", "H3N2", "H3N2",
                 "Influenza B", "Influenza B", "Influenza B")

virusname = c("h1n1sea_sub1","h1n1sea_sub2","h1n1sea_sub3",
              "h1n1pdm_sub1","h1n1pdm_sub2","h1n1pdm_sub3", 
              "h3n2_sub1", "h3n2_sub2", "h3n2_sub3", 
              "infB_sub1", "infB_sub2", "infB_sub3")

replicate = c("1","2","3",
              "1","2","3", 
              "1", "2", "3", 
              "1", "2", "3")

tick_labels = c("p1918 like H1N1",
                "p2009 like H1N1",
                "H3N2",
                "Influenza B")

offset=0.25

xval = c(1-offset,1,1+offset,
         2-offset,2,2+offset,
         3-offset,3,3+offset,
         4-offset,4,4+offset)


# read in the clade support and clade height
first = T;
for (i in seq(1, length(viruses),3)){
  fname = paste("./", viruses[[i]], "/trees/Heights.csv", sep="")
  t = read.table(fname, header=T, sep=",", na.strings="")
  
  fname = paste("./", viruses[[i]], "/trees/Clades.csv", sep="")
  t.Clades = read.table(fname, header=T, sep=",", na.strings="")
  
  t = t[-which(t$post1<0.5 | t$post2<0.5),]
  t$virtype = viruses[[i]]

  if (first){
    heights = t
    clades = t.Clades
    first=F
  }else{
    heights = rbind(heights, t)
    clades = rbind(clades, t.Clades)
  }
}

heights$ratio = ((heights$upper2-heights$lower2)/heights$median2)/ ((heights$upper1-heights$lower1)/heights$median1)
heights$x = heights$ratio
for (i in seq(1,length(virusname))){
  heights[which(heights$virus==virusname[[i]]), "x"] = xval[[i]]
  clades[which(heights$virus==virusname[[i]]), "x"] = xval[[i]]
}
for (i in seq(1, length(viruses),3)){
  med.vals = median(heights[which(heights$virtype== viruses[[i]]), "ratio"])
  
  if (i==1){
    median.vals = data.frame(x = xval[[i+1]], y = 3, text=format(med.vals,digits=2))
  }else{
    median.vals = rbind(median.vals, data.frame(x = xval[[i+1]], y = 3, text=format(med.vals,digits=2)))
  }
}
# get the medians

p.heights <- ggplot(heights) +
  geom_violin(aes(x=x, y=ratio, group=x),color="black", fill="white") +
  geom_text(data=median.vals, aes(x=x,y=y,label=text))+
  scale_y_log10(limits=c(0.1,10), breaks=c(0.1,0.5, 1.0,2,10)) +
  theme_minimal()+
  ylab("ratio of relative node height uncertainty") +
  xlab("")+
  scale_x_continuous(breaks=seq(1,4), label=tick_labels) 
plot(p.heights)

clades$label = clades$method

p.clades <- ggplot(clades) +
  geom_violin(aes(x=x, y=post, group=interaction(x,method), fill=method, color=method)) +
  scale_y_continuous(limits=c(0,1)) +
  theme_minimal()+
  ylab("posterior node support") +
  xlab("")+
  scale_fill_brewer(type='qual') +
  scale_color_brewer(type='qual') +
  scale_x_continuous(breaks=seq(1,4), label=tick_labels) +
  theme(legend.position = "none")
plot(p.clades)



first = T;
for(i in seq(1, length(viruses))){
  # set the log file name
  filename1 = path=paste("./",viruses[[i]] ,"/combined/", virusname[[i]],".log", sep="")
  filename2 = path=paste("./",viruses[[i]] ,"/combined/", virusname[[i]], "norea",".log", sep="")
  
  # read the file in
  t <- read.table(filename1, header=TRUE, sep="\t")
  t.norea <- read.table(filename2, header=TRUE, sep="\t")
  
  # get the HPD's
  hpd.ne = HPDinterval(as.mcmc(t$popSize.t))
  hpd.clock = HPDinterval(as.mcmc(t$clockRate.c))
  
  hpd.ne.norea = HPDinterval(as.mcmc(t.norea$popSize.t))
  hpd.clock.norea = HPDinterval(as.mcmc(t.norea$clockRate.c))
  
  
  new.reassortment = data.frame(ne.lower=hpd.ne[1,"lower"], ne.upper=hpd.ne[1,"upper"],
                                clock.lower=hpd.clock[1,"lower"], clock.upper=hpd.clock[1,"upper"])
  
  new.reassortment.norea = data.frame(ne.lower=hpd.ne.norea[1,"lower"], ne.upper=hpd.ne.norea[1,"upper"],
                                clock.lower=hpd.clock.norea[1,"lower"], clock.upper=hpd.clock.norea[1,"upper"])
  
  new.reassortment$virus = virus_labels[[i]]
  new.reassortment$replicate = replicate[[i]]
  new.reassortment$xval = xval[[i]]
  
  new.reassortment.norea$virus = virus_labels[[i]]
  new.reassortment.norea$replicate = replicate[[i]]
  new.reassortment.norea$xval = xval[[i]]
  
  
  if (first){
    reassortment = new.reassortment
    reassortment.norea = new.reassortment.norea
    first = F;
  }else{
    reassortment = rbind(reassortment, new.reassortment)
    reassortment.norea = rbind(reassortment.norea, new.reassortment.norea)
  }
}





p.ne <- ggplot() +
  geom_linerange(data=reassortment, aes(x=xval, ymin=ne.lower, ymax=ne.upper, color=" with reassortment"), size=4) +
  geom_linerange(data=reassortment.norea, aes(x=xval, ymin=ne.lower, ymax=ne.upper, color=" independent segments"), size=4) +
  xlab("") +
  ylab("effective population size") +
  theme_minimal() +
  scale_x_continuous(breaks=seq(1,4), label=tick_labels) +
  scale_color_brewer(type='qual') +
  theme(legend.position=c(0.4, 0.8),
        legend.title = element_blank(),
        legend.text = element_text(size=15))
plot(p.ne)


p.clock <- ggplot(reassortment) +
  geom_linerange(data=reassortment, aes(x=xval, ymin=clock.lower, ymax=clock.upper, color=" with reassortment"), size=4) +
  geom_linerange(data=reassortment.norea, aes(x=xval, ymin=clock.lower, ymax=clock.upper, color=" independent"), size=4) +
  xlab("") +
  ylab("evolutionary rate pre site and year") +
  theme_minimal() +
  scale_x_continuous(breaks=seq(1,4), label=tick_labels) +
  scale_color_brewer(type='qual') +
  theme(legend.position="none")
        
plot(p.clock)

plot.all = grid.arrange(p.heights, p.clades, p.ne, p.clock, ncol=2)


ggsave(plot=plot.all,paste("../../Reassortment-Text/Figures/precision_bias.pdf", sep=""),width=10, height=6)

