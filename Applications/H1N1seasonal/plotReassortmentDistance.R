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
threshold = 1

# define the number of segments
nrSegments = 9

time =""

# read in the file with the reassortment distances
fileName= paste("combined/h1n1sea", time, ".distance.txt", sep="")
con=file(fileName,open="r")
line=readLines(con) 
first = T
for (i in seq(1,length(line))){
  print(paste(i, "of", length(line)) )
  
  splitline = strsplit(line[[i]], split="\t")
  firstevent=T
  remove(dat)
  c=1
  
  for (j in seq(2,length(splitline[[1]]))){
    tmp = strsplit(splitline[[1]][[j]], split=",")
    tmp2 = strsplit(tmp[[1]][[1]], split=":")
    # checks if the event is on the trunk
    if (as.logical(tmp2[[1]][[1]])){
      c=c+1
      for (k in seq(2,length(tmp[[1]]))){
        tmp2 = strsplit(tmp[[1]][[k]], split=":")
          if (!is.na(as.numeric(tmp2[[1]][[3]]))){
            distance = as.numeric(tmp2[[1]][[3]])-as.numeric(tmp2[[1]][[2]])
            if (distance>threshold){
            tmp3 = strsplit(tmp2[[1]][[1]], split="-")
            seg1 = as.numeric(tmp3[[1]][[1]])
            seg2 = as.numeric(tmp3[[1]][[2]])
            new.dat = data.frame(from=seg1,to=seg2)
            if(firstevent){
              dat = new.dat
              firstevent = F
            }else{
              dat = rbind(dat,new.dat)
            }
          }
        }
      }
    # }
  }
  
  # count how many events there were between the segments
  for (a in seq(0,nrSegments-2)){
    for (b in seq(a+1,nrSegments-2)){
      if (a!=b){
        obs = which(dat$from==a & dat$to==b)
        obs = union(obs, which(dat$from==b & dat$to==a))
        
        prior = which(dat$from==a & dat$to==nrSegments-1)
        prior = union(prior, which(dat$from==b & dat$to==nrSegments-1))

        prior = union(prior, which(dat$from==nrSegments-1 & dat$to==a))
        prior = union(prior, which(dat$from==nrSegments-1 & dat$to==b))
        
        
        
        # ks.results.less = ks.test(dat[obs,"dist"], dat[prior,"dist"], alternative = "less")
        # ks.results.more = ks.test(dat[obs,"dist"], dat[prior,"dist"], alternative = "greater")
        if (2*length(obs)>length(prior)){
          new.reascount = data.frame(from=a,to=b,larger=1, smaller=0, obs=length(obs))
        }else if (2*length(obs)<length(prior)){
          new.reascount = data.frame(from=a,to=b,larger=0, smaller=1, obs=length(obs))
        }else{
          new.reascount = data.frame(from=a,to=b,larger=0, smaller=0, obs=length(obs))
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
  
  # # count how many events there were between the segments
  # for (a in seq(0,nrSegments-2)){
  #   for (b in seq(0,nrSegments-2)){
  #     if (a!=b){
  #       obs = which(dat$from==a & dat$to==b)
  #       prior = which(dat$from==a & dat$to==nrSegments-1)
  #       ks.results.less = ks.test(dat[obs,"dist"], dat[prior,"dist"], alternative = "less")
  #       ks.results.more = ks.test(dat[obs,"dist"], dat[prior,"dist"], alternative = "greater")
  #       if (ks.results.less$p.value<ks.results.more$p.value){
  #         new.reascount = data.frame(from=a,to=b,larger=1, sum.diff = sum(dat[obs,"dist"])-sum(dat[prior,"dist"]),
  #                                    max.dist = max(dat[obs,"dist"]), max.dist.prior = max(dat[prior,"dist"]) )
  #       }else{
  #         new.reascount = data.frame(from=a,to=b,larger=0, sum.diff = sum(dat[obs,"dist"])-sum(dat[prior,"dist"]),
  #                                    max.dist = max(dat[obs,"dist"]), max.dist.prior = max(dat[prior,"dist"]))
  #       }
  #       if(first){
  #         reascount = new.reascount
  #         first = F
  #       }else{
  #         reascount = rbind(reascount, new.reascount)
  #       }
  #     }
  #   }
  # }
  
  
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

# reascount$from = segments[reascount$from+1]
# reascount$to = segments[reascount$to+1]

# plot the histograms
p_ks = ggplot()+
  geom_tile(data=dat, aes(x=from, y=to, fill=mean))+
  # facet_grid(.~to) +
  # scale_fill_OkabeIto() +
  scale_fill_gradientn(colors = c("#88419d","#b3cde3",
                                  rgb(1,1,1),  
                                  "#fdcc8a","#d7301f"),
                       limits = c(0,1), 
                       na.value = "white",
                       breaks = c(0, 0.1, 0.5, 0.9, 1), 
                       name="obs>prior") +
  xlab("") +
  ylab("") +
  theme_minimal()
plot(p_ks)


ggsave(plot=p_ks,paste("../../../Reassortment-Text/Figures/distance/H1N1sea", time, ".pdf", sep=""), width=5, height=4)


indices = which(reascount$from==1 & reascount$to==2)
reascount[indices, "obs"]