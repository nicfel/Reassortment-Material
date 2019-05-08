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

# define the threshold of when to count an event
mrca.threshold = 35

# define the number of segments
nrSegments = 9

mrsi = 2017

# read in the file with the reassortment distances
fileName="combined/infB.distance.txt"
con=file(fileName,open="r")
line=readLines(con) 
print("loaded file")

first = T
for (i in seq(1,length(line))){
  print(i)

  splitline = strsplit(line[[i]], split="\t")
  firstevent=T
  remove(dat)
  for (j in seq(2,length(splitline[[1]]))){
    tmp = strsplit(splitline[[1]][[j]], split=",")
    for (k in seq(1,length(tmp[[1]]))){
      
      tmp2 = strsplit(tmp[[1]][[k]], split=":")
      

      tmp3 = strsplit(tmp2[[1]][[1]], split="-")
      
      seg1 = as.numeric(tmp3[[1]][[1]])
      seg2 = as.numeric(tmp3[[1]][[2]])
      
      new.dat = data.frame(from=min(seg1,seg2),to=max(seg1,seg2), mrca = as.numeric(tmp2[[1]][[3]]), nodeheight=as.numeric(tmp2[[1]][[2]]))
      
      if(firstevent){
        dat = new.dat
        firstevent = F
      }else{
        dat = rbind(dat,new.dat)
      }
    }
  }

  # get the minimal node height of an event that for which the mrca is larget than a threshold
  for (a in seq(0,0)){
    for (b in seq(a+1,nrSegments-1)){
    # for (b in seq(2,2)){
        
      indices = which(dat$from==a & dat$to==b & dat$mrca>mrca.threshold)
      # for (j in seq(1,length(indices))){

        new.reascount = data.frame(from=a,to=b, min=min(dat[indices,"nodeheight"]))
        
        if(first){
          reascount = new.reascount
          first = F
        }else{
          reascount = rbind(reascount, new.reascount)
        }
      }
    # }
  }
}



segments = c("HA", "MP", "NA", "NP", "NS", "PA", "PB1", "PB2", "prior")
reascount$from = segments[reascount$from+1]
reascount$to = segments[reascount$to+1]


# plot the histograms

p_hist = ggplot()+
  geom_violin(data=reascount, aes(x=to, y=2016-min))+
  theme_minimal()
plot(p_hist)

