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

# read in the file with the reassortment distances
fileName="combined/h1n1sea.distance.txt"
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
        obs = length(which(dat[which(dat$from==a),]$to==b))
        prior = length(which(dat[which(dat$from==a),]$to==nrSegments-1))
        new.reascount = data.frame(from=a,to=b,obs=obs, prior=prior)
        # new.reascount = rbind(new.reascount, data.frame(from=b,to=a,obs=obs, prior=prior))
        new.reascount.prior = data.frame(from=a,to=b,obs=prior)
        # new.reascount.prior = rbind(new.reascount.prior, data.frame(from=b,to=a,obs=prior))
        if(first){
          reascount = new.reascount
          reascount.prior = new.reascount.prior
          first = F
        }else{
          reascount = rbind(reascount, new.reascount)
          reascount.prior = rbind(reascount.prior, new.reascount.prior)
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
      difference = reascount[indices,"obs"] - reascount[indices,"prior"]
      hpd = HPDinterval(as.mcmc(difference))
      new.dat = data.frame(from=a,to=b,mean=mean(difference), median=median(difference),upper=hpd[1,"upper"],lower=hpd[1,"lower"])
      if (first){
        dat=new.dat
        first = F
      }else{
        dat=rbind(dat,new.dat)
      }
    }
  }
}

segments = c("HA", "MP", "NA", "NP", "NS", "PA", "PB1", "PB2", "prior")
dat$from = segments[dat$from+1]
dat$to = segments[dat$to+1]




# plot the histograms

p_hist = ggplot()+
  geom_crossbar(data=dat, aes(x=from, y=mean,ymin=lower, ymax=upper,fill=from))+
  geom_hline(yintercept=0, linetype="dashed", 
             color = "red", size=2)+
  facet_grid(.~to) +
  scale_fill_OkabeIto() +
  theme_light()
plot(p_hist)

