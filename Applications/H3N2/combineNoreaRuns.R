######################################################
######################################################
# combine the gene tree runs and run the mcc trees
######################################################
######################################################
# clear workspace
rm(list = ls())

# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)


trees <- list.files(path=paste("./out/",sep=""), pattern="*norea_rep0.*\\.log$", full.names = TRUE)

for (i in seq(1,length(trees))){
  in_command <- " -b 20 -resample 100000 -log"
  for (j in seq(0,2)){
    in_command = paste(in_command, " ", gsub("rep0", paste("rep", j,sep=""), trees[i]), sep="")
  }
  
  out_command = gsub("rep0_", "", trees[i])
  out_command = gsub("out", "combined", out_command)
  
  combined_command = gsub(".trees",".trees", out_command)
  combined_command = paste(" -o ", gsub("_rep0", "",combined_command), sep="")
  # combine the trees
  system(paste("/Applications/BEAST\\ 2.5.2/bin/logcombiner", gsub(".network.trees",".log", in_command), gsub(".network.trees",".log", combined_command), sep=" "))
}

trees <- list.files(path=paste("./out/",sep=""), pattern="*ind_rep0.*\\.log$", full.names = TRUE)

for (i in seq(1,length(trees))){
  in_command <- " -b 20 -resample 100000 -log"
  for (j in seq(0,2)){
    in_command = paste(in_command, " ", gsub("rep0", paste("rep", j,sep=""), trees[i]), sep="")
  }
  
  out_command = gsub("rep0_", "", trees[i])
  out_command = gsub("out", "combined", out_command)
  
  combined_command = gsub(".trees",".trees", out_command)
  combined_command = paste(" -o ", gsub("_rep0", "",combined_command), sep="")
  # combine the trees
  system(paste("/Applications/BEAST\\ 2.5.2/bin/logcombiner", gsub(".network.trees",".log", in_command), gsub(".network.trees",".log", combined_command), sep=" "))
}


