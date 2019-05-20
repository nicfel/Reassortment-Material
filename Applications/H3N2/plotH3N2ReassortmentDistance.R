######################################################
######################################################
# Here the inferred mean coalescent and migration
# rate ratios are plotted
######################################################
######################################################
library(ggplot2)
library(coda)
library(colorblindr)

# clear workspace
rm(list = ls())

# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

# define the threshold of when to count an event
threshold = 0

# define the number of segments
nrSegments = 9

time =""

# read in the file with the reassortment distances
fileName= paste("combined/h3n2", time, ".distance.txt", sep="")
con=file(fileName,open="r")
line=readLines(con) 
first = T
for (i in seq(1000,length(line))){
  print(paste(i, "of", length(line)) )
  
  splitline = strsplit(line[[i]], split="\t")
  # firstevent=T
  remove(dat)
  c=1
  
  firstevent = T
  for (j in seq(2,length(splitline[[1]]))){
    tmp = strsplit(splitline[[1]][[j]], split=",")
    tmp2 = strsplit(tmp[[1]][[1]], split=":")
    # checks if the event is on the trunk
    isTrunk = as.logical(tmp2[[1]][[1]])
    c=c+1
    firstseg = T
    for (k in seq(2,length(tmp[[1]]))){
      tmp2 = strsplit(tmp[[1]][[k]], split=":")
      if (!is.na(as.numeric(tmp2[[1]][[3]]))){
        distance = as.numeric(tmp2[[1]][[3]])-as.numeric(tmp2[[1]][[2]])
        if (distance>threshold){
          tmp3 = strsplit(tmp2[[1]][[1]], split="-")
          seg1 = as.numeric(tmp3[[1]][[1]])
          seg2 = as.numeric(tmp3[[1]][[2]])
          new.seg.dat = data.frame(from=seg1,to=seg2, dist = distance, trunk = isTrunk )
          
          if(firstevent){
            seg.dat = new.seg.dat
            firstevent = F
          }else{
            seg.dat = rbind(seg.dat,new.seg.dat)
          }# firstseg
        }# if threshold
      }# if non na
      #segments
    }#events
  }
  # compute the ecdf
  for (a in seq(0,nrSegments-1)){
    for (b in seq(0,nrSegments-1)){
      if (a!=b){
        obs.trunk = which(seg.dat$trunk & seg.dat$from==a & seg.dat$to==b)
        obs.off.trunk = which(!seg.dat$trunk & seg.dat$from==a & seg.dat$to==b)
        e.trunk = ecdf(seg.dat[obs.trunk,"dist"])
        e.trunk.off = ecdf(seg.dat[obs.off.trunk,"dist"])
        
        new.ecdf = data.frame(seg1=a,seg2=b,
                              x=seg.dat[obs.trunk,"dist"],
                              y=e.trunk(seg.dat[obs.trunk,"dist"]), trunk="is trunk")
        new.ecdf = rbind(new.ecdf, data.frame(seg1=a,seg2=b,
                                              x=seg.dat[obs.off.trunk,"dist"],
                                              y=e.trunk.off(seg.dat[obs.off.trunk,"dist"]), trunk="off trunk"))
        if(first){
          # reascount = new.reascount
          ecdf.dat = new.ecdf
          first = F
        }else{
          # reascount = rbind(reascount, new.reascount)
          ecdf.dat = rbind(ecdf.dat, new.ecdf)
        }
        
      }
    }
  }
  
  
  # ks.results.less = ks.test(dat[obs.trunk,"dist"], dat[obs.off.trunk,"dist"], alternative = "less")
  # ks.results.more = ks.test(dat[obs.trunk,"dist"], dat[obs.off.trunk,"dist"], alternative = "greater")
  # 
  # if (ks.results.more$p.value<ks.results.less$p.value){
  #   new.reascount = data.frame(larger=1)
  # }else{
  #   new.reascount = data.frame(larger=0)
  # }
  # 
  
}

segments = c("HA", "MP", "NA", "NP", "NS", "PA", "PB1", "PB2", "dummy")
ecdf.dat$seg1 = segments[ecdf.dat$seg1+1]
ecdf.dat$seg2 = segments[ecdf.dat$seg2+1]

ecdf.dat$seg1 <- factor(ecdf.dat$seg1, levels =segments)
ecdf.dat$seg2 <- factor(ecdf.dat$seg2, levels =segments)


p <- ggplot(ecdf.dat) + geom_point(aes(x=x,y=y,color=trunk), alpha=0.05)+
  geom_smooth(aes(x=x,y=y,color=trunk)) +
  scale_x_continuous(limits=c(0,7.5)) +
  scale_y_continuous(limits=c(0,1)) +
  facet_grid(seg1~seg2)+
  xlab("distances in years")+
  ylab("cumulative probability (empirical)")+
  theme_light()+
  scale_color_OkabeIto()
plot(p)
ggsave(plot=p,paste("../../../Reassortment-Text/Figures/distance/H3N2.distance.pdf", sep=""), width=10, height=10)

