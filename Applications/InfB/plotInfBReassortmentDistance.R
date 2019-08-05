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
fileName= paste("combined/infB_sub2", time, ".distance.txt", sep="")
con=file(fileName,open="r")
line=readLines(con) 
first = T
for (i in seq(800,length(line))){
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
    }#events
  }
  
  # compute the ecdf
  for (a in seq(0,nrSegments-2)){
    for (b in seq(0,nrSegments-2)){
      if (a!=b){
        obs.trunk = which(seg.dat$trunk & seg.dat$from==a & seg.dat$to==b)
        obs.off.trunk = which(!seg.dat$trunk & seg.dat$from==a & seg.dat$to==b)
        
        e.trunk = seg.dat[obs.trunk,"dist"]
        e.trunk.off = seg.dat[obs.off.trunk,"dist"]
        
        obs.trunk.dummy = which(seg.dat$trunk & seg.dat$from==a & seg.dat$to==nrSegments-1)
        obs.off.trunk.dummy = which(!seg.dat$trunk & seg.dat$from==a & seg.dat$to==nrSegments-1)
        
        e.trunk.dummy = seg.dat[obs.trunk.dummy,"dist"]
        e.trunk.off.dummy = seg.dat[obs.off.trunk.dummy,"dist"]
        
        
        # new.ecdf = data.frame(seg1=a,seg2=b,
        #                       x=seg.dat[obs.trunk,"dist"],
        #                       y=e.trunk(seg.dat[obs.trunk,"dist"]), trunk="is trunk")
        # new.ecdf = rbind(new.ecdf, data.frame(seg1=a,seg2=b,
        #                                       x=seg.dat[obs.off.trunk,"dist"],
        #                                       y=e.trunk.off(seg.dat[obs.off.trunk,"dist"]), trunk="off trunk"))
        
        new.ecdf = data.frame(seg1=a,seg2=b,
                              nr.trunk=length(e.trunk), nr.nontrunk=length(e.trunk.off),
                              nr.trunk.dummy=length(e.trunk.dummy), nr.nontrunk.dummy=length(e.trunk.off.dummy))
        
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
}
segments = c("HA", "MP", "NA", "NP", "NS", "PA", "PB1", "PB2", "dummy")
ecdf.dat$from = segments[ecdf.dat$seg1+1]
ecdf.dat$to = segments[ecdf.dat$seg2+1]

ecdf.dat$from <- factor(ecdf.dat$from, levels =segments)
ecdf.dat$to <- factor(ecdf.dat$to, levels =segments)

p <-ggplot(ecdf.dat) + 
  geom_histogram(aes( (nr.trunk+nr.nontrunk) -  (nr.trunk.dummy+nr.nontrunk.dummy)), color=NA, fill="black", binwidth=1)+
  geom_vline(xintercept = 0, color="red")+
  facet_grid(from~to)+
  theme_light()+
  scale_color_OkabeIto()
plot(p)
ggsave(plot=p,paste("../../../Reassortment-Text/Figures/distance/InfB_sub2_nrevents.pdf", sep=""), width=10, height=10)
dsa

p.trunk <- ggplot(ecdf.dat) + 
  geom_histogram(aes( (nr.trunk+nr.nontrunk) -  (nr.trunk.dummy+nr.nontrunk.dummy)), fill="black", binwidth=1, alpha=0.2)+
  geom_histogram(aes( (nr.trunk) -  (nr.trunk.dummy), fill="fit"), binwidth=1, alpha=0.2)+
  geom_histogram(aes( (nr.nontrunk) -  (nr.nontrunk.dummy), fill="unfit"), binwidth=1, alpha=0.2)+
  geom_vline(xintercept = 0)+
  facet_grid(from~to)+
  theme_light()+
  scale_fill_OkabeIto() +
  ggtitle("fit events")
plot(p.trunk)

p.offtrunk <-ggplot(ecdf.dat) + 
  geom_histogram(aes( (nr.nontrunk) -  (nr.nontrunk.dummy)), binwidth=1)+
  geom_vline(xintercept = 0)+
  facet_grid(from~to)+
  theme_light()+
  scale_color_OkabeIto() +
  ggtitle("unfit events")
plot(p.offtrunk)




# p <- ggplot(ecdf.dat) + 
#   geom_violin(aes(x=1,y=p)) +
#   # scale_x_continuous(limits=c(0,7.5)) +
#   scale_y_log10(limits=c(0.001,1)) +
#   facet_grid(from~to)+
#   xlab("distances in years")+
#   ylab("cumulative probability (empirical)")+
#   theme_light()+
#   scale_color_OkabeIto()
# plot(p)
# # ggsave(plot=p,paste("../../../Reassortment-Text/Figures/distance/H3N2_distance.pdf", sep=""), width=10, height=10)

