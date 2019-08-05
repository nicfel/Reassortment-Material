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

fnames = c("h3n2recent")
for (fn in seq(1,length(fnames))){
  # read in the file with the reassortment distances
  fileName= paste("combined/", fnames[[fn]], ".distance.txt", sep="")
  con=file(fileName,open="r")
  line=readLines(con) 
  first = T
  for (i in seq(1091,length(line))){
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
  
  first = T
  # build heatmap
  for (a in seq(0,nrSegments-2)){
    for (b in seq(0,nrSegments-2)){
      if (a!=b){
        indices = which(ecdf.dat$seg1==a & ecdf.dat$seg2==b)

        all = (ecdf.dat[indices, "nr.trunk"] + ecdf.dat[indices, "nr.nontrunk"]) - (ecdf.dat[indices, "nr.trunk.dummy"] + ecdf.dat[indices, "nr.nontrunk.dummy"])
        fit = ecdf.dat[indices, "nr.trunk"]- ecdf.dat[indices, "nr.trunk.dummy"]
        unfit = ecdf.dat[indices, "nr.nontrunk"] - ecdf.dat[indices, "nr.nontrunk.dummy"]
        
        if (median(all)>0){new.plot.vals.all = data.frame(seg1=a,seg2=b,val=-log10(sum(all<=0)/length(all)))
        }else{new.plot.vals.all = data.frame(seg1=a,seg2=b,val=log10(sum(all>=0)/length(all)))}

        if (median(all)>0){new.plot.vals.unfit = data.frame(seg1=a,seg2=b,val=1-sum(fit<=0)/length(fit))
        }else{new.plot.vals.unfit = data.frame(seg1=a,seg2=b,val=-1+sum(fit>=0)/length(fit))}
        
        if (median(all)>0){new.plot.vals.fit = data.frame(seg1=a,seg2=b,val=1-sum(unfit<=0)/length(unfit))
        }else{new.plot.vals.fit = data.frame(seg1=a,seg2=b,val=-1+sum(unfit>=0)/length(unfit))}
        
      }else{
        new.plot.vals.all = data.frame(seg1=a,seg2=b,val=NA)
        new.plot.vals.unfit = data.frame(seg1=a,seg2=b,val=NA)
        new.plot.vals.fit = data.frame(seg1=a,seg2=b,val=NA)
      }

      if(first){
        plot.vals.all = new.plot.vals.all
        plot.vals.fit = new.plot.vals.unfit
        plot.vals.unfit = new.plot.vals.fit
        first = F
      }else{
        plot.vals.all = rbind(new.plot.vals.all, plot.vals.all)
        plot.vals.fit = rbind(new.plot.vals.unfit, plot.vals.fit)
        plot.vals.unfit = rbind(new.plot.vals.fit, plot.vals.unfit)
      }
      
    }
  }
        
  segments = c("HA", "MP", "NA", "NP", "NS", "PA", "PB1", "PB2", "dummy")
  
  plot.vals.all$from = segments[plot.vals.all$seg1+1]
  plot.vals.all$to = segments[plot.vals.all$seg2+1]
  plot.vals.all$from <- factor(plot.vals.all$from, levels =segments)
  plot.vals.all$to <- factor(plot.vals.all$to, levels =segments)
  
  plot.vals.fit$from = segments[plot.vals.fit$seg1+1]
  plot.vals.fit$to = segments[plot.vals.fit$seg2+1]
  plot.vals.fit$from <- factor(plot.vals.fit$from, levels =segments)
  plot.vals.fit$to <- factor(plot.vals.fit$to, levels =segments)
  
  plot.vals.unfit$from = segments[plot.vals.unfit$seg1+1]
  plot.vals.unfit$to = segments[plot.vals.unfit$seg2+1]
  plot.vals.unfit$from <- factor(plot.vals.unfit$from, levels =segments)
  plot.vals.unfit$to <- factor(plot.vals.unfit$to, levels =segments)
  
  
  
  
  # 
  # 
  # 
  # 
  # p <-ggplot(ecdf.dat) + 
  #   geom_histogram(aes( (nr.trunk+nr.nontrunk) -  (nr.trunk.dummy+nr.nontrunk.dummy)), binwidth=1)+
  #   geom_vline(xintercept = 0)+
  #   facet_grid(from~to)+
  #   theme_light()+
  #   scale_color_OkabeIto()
  # plot(p)
  
  lim=2
  plot.vals.all$val[which(plot.vals.all$val< -lim)] <- -lim
  plot.vals.all$val[which(plot.vals.all$val> lim)] <- lim
  
  p <- ggplot(plot.vals.all)+
    geom_tile(aes(x=from,y=to,fill=val))+
    theme_minimal()+
    scale_fill_gradientn(na.value="white", colors = c("#88419d","#b3cde3", rgb(1,1,1),  "#fdcc8a","#d7301f"),
                         limits = c(-lim,lim), 
                         breaks = c(-1, -0.5, 0, 0.5, 1), 
                         name="Association",
                         label=c(expression(paste("p" <= 0.01, sep="")),"p=0.1","","p=0.1",expression(paste("p" <= 0.01, sep=""))
                         ))
  plot(p)
    
  ggsave(plot=p,paste("../../../Reassortment-Text/Figures/distance/", fnames[[fn]], "_nrevents.pdf", sep=""), width=3, height=3)
  
  
  p <- ggplot(plot.vals.fit)+
    geom_tile(aes(x=from,y=to,fill=val))+
    theme_minimal()+
    scale_fill_gradientn(na.value="white", colors = c("#88419d","#b3cde3", rgb(1,1,1),  "#fdcc8a","#d7301f"),
                         limits = c(-lim,lim), 
                         breaks = c(-1, -0.5, 0, 0.5, 1), 
                         name="Association",
                         label=c(expression(paste("p" <= 0.01, sep="")),"p=0.1","","p=0.1",expression(paste("p" <= 0.01, sep=""))
                         ))
  plot(p)
  ggsave(plot=p,paste("../../../Reassortment-Text/Figures/distance/", fnames[[fn]], "_nrevents.fit.pdf", sep=""), width=3, height=3)
  
  
  p <- ggplot(plot.vals.unfit)+
    geom_tile(aes(x=from,y=to,fill=val))+
    theme_minimal()+
    scale_fill_gradientn(na.value="white", colors = c("#88419d","#b3cde3", rgb(1,1,1),  "#fdcc8a","#d7301f"),
                         limits = c(-lim,lim), 
                         breaks = c(-1, -0.5, 0, 0.5, 1), 
                         name="Association",
                         label=c(expression(paste("p" <= 0.01, sep="")),"p=0.1","","p=0.1",expression(paste("p" <= 0.01, sep=""))
                         ))
  plot(p)
  ggsave(plot=p,paste("../../../Reassortment-Text/Figures/distance/", fnames[[fn]], "_nrevents.unfit.pdf", sep=""), width=3, height=3)
  
  
  # 
  # p.trunk <-ggplot(ecdf.dat) + 
  #   geom_histogram(aes( (nr.trunk) -  (nr.trunk.dummy)), binwidth=1)+
  #   geom_vline(xintercept = 0)+
  #   facet_grid(from~to)+
  #   theme_light()+
  #   scale_color_OkabeIto()
  # plot(p.trunk)
  # ggsave(plot=p.trunk,paste("../../../Reassortment-Text/Figures/distance/", fnames[[fn]], "_nrevents.fit.pdf", sep=""), width=10, height=10)
  # 
  # p.offtrunk <-ggplot(ecdf.dat) + 
  #   geom_histogram(aes( (nr.nontrunk) -  (nr.nontrunk.dummy)), binwidth=1)+
  #   geom_vline(xintercept = 0)+
  #   facet_grid(from~to)+
  #   theme_light()+
  #   scale_color_OkabeIto()
  # plot(p.offtrunk)
  # ggsave(plot=p.offtrunk,paste("../../../Reassortment-Text/Figures/distance/", fnames[[fn]], "_nrevents.unfit.pdf", sep=""), width=10, height=10)

}


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

