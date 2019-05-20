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

# read in the file with the reassortment distances
fileName="combined/h2n2.distance.txt"
con=file(fileName,open="r")
line=readLines(con) 
first = T
for (i in seq(1,length(line))){
  print(paste(i, "of", length(line)) )
  
  splitline = strsplit(line[[i]], split="\t")
  firstevent=T
  remove(dat)
  
  for (j in seq(2,length(splitline[[1]]))){
    tmp = strsplit(splitline[[1]][[j]], split=",")
    for (k in seq(1,length(tmp[[1]]))){
      tmp2 = strsplit(tmp[[1]][[k]], split=":")
      if (!is.na(as.numeric(tmp2[[1]][[3]]))){
        distance = as.numeric(tmp2[[1]][[3]])-as.numeric(tmp2[[1]][[2]])
        if (distance>threshold){
          tmp3 = strsplit(tmp2[[1]][[1]], split="-")
          seg1 = as.numeric(tmp3[[1]][[1]])
          seg2 = as.numeric(tmp3[[1]][[2]])
          new.dat = data.frame(from=seg1,to=seg2, dist=distance)
          if(firstevent){
            dat = new.dat
            firstevent = F
          }else{
            dat = rbind(dat,new.dat)
          }
        }
      }
    }
  }
  # count how many events there were between the segments
  for (a in seq(0,nrSegments-2)){
    for (b in seq(0,nrSegments-2)){
      if (a!=b){
        obs = which(dat$from==a & dat$to==b)
        prior = which(dat$from==a & dat$to==nrSegments-1)
        ks.results.less = ks.test(dat[obs,"dist"], dat[prior,"dist"], alternative = "less")
        ks.results.more = ks.test(dat[obs,"dist"], dat[prior,"dist"], alternative = "greater")
        if (ks.results.less$p.value<ks.results.more$p.value){
            new.reascount = data.frame(from=a,to=b,larger=1, sum.diff = sum(dat[obs,"dist"])-sum(dat[prior,"dist"]),
                                   max.dist = max(dat[obs,"dist"]), max.dist.prior = max(dat[prior,"dist"]) )
        }else{
            new.reascount = data.frame(from=a,to=b,larger=0, sum.diff = sum(dat[obs,"dist"])-sum(dat[prior,"dist"]),
                                     max.dist = max(dat[obs,"dist"]), max.dist.prior = max(dat[prior,"dist"]))
        }
        if(first){
          reascount = new.reascount
          first = F
        }else{
          reascount = rbind(reascount, new.reascount)
        }
      }
    }
  }
  
}

# get the 95% HPD
first = T
for (a in seq(0,nrSegments-2)){
  for (b in seq(0,nrSegments-2)){
    if (a!=b){
      indices = which(reascount$from==a & reascount$to==b)
      indices = union(indices, which(reascount$from==b & reascount$to==a))
      new.dat = data.frame(from=a,to=b,mean=mean(reascount[indices,"larger"]))
      if (first){
        dat=new.dat
        # dat.p=new.dat.p
        first = F
      }else{
        dat=rbind(dat,new.dat)
        # dat.p=rbind(dat.p,new.dat.p)
      }
    }else{
      new.dat = data.frame(from=a,to=b,mean=NA)
      if (first){
        dat=new.dat
        # dat.p=new.dat.p
        first = F
      }else{
        dat=rbind(dat,new.dat)
        # dat.p=rbind(dat.p,new.dat.p)
      }
    }
  }
}


segments = c("HA", "MP", "NA", "NP", "NS", "PA", "PB1", "PB2", "prior")
dat$from = segments[dat$from+1]
dat$to = segments[dat$to+1]

reascount$from = segments[reascount$from+1]
reascount$to = segments[reascount$to+1]



# plot the histograms

p_ks = ggplot()+
  geom_tile(data=dat, aes(x=from, y=to, fill=mean))+
  # facet_grid(.~to) +
  # scale_fill_OkabeIto() +
  scale_fill_gradientn(values=seq(0,1,0.25), breaks=seq(0,1,0.25),  colours=c("black", "white", "white", "white", "red") , na.value="white")+
  xlab("") +
  ylab("") +
  theme_minimal()
plot(p_ks)


ggsave(plot=p_ks,paste("../../../Reassortment-Text/Figures/distance/H2N2.pdf", sep=""), width=5, height=4)



# p_ks = ggplot()+
#   geom_violin(data=reascount, aes(x=from, y=ks, color=from))+
#   facet_grid(from~to) +
#   scale_fill_OkabeIto() +
#   geom_hline(yintercept=0, 
#              color = "red", size=1)+
#   theme_light()
# plot(p_ks)


# p_p = ggplot()+
#   geom_crossbar(data=dat.p, aes(x=from, y=mean,ymin=lower, ymax=upper,fill=from))+
#   facet_grid(.~to) +
#   scale_fill_OkabeIto() +
#   geom_hline(yintercept=0, 
#              color = "red", size=1)+
#   theme_light()
# plot(p_p)


