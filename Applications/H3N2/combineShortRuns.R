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


system("rm -r shortout/combined")
system("mkdir shortout/combined")


trees <- list.files(path=paste("./shortout/individual/",sep=""), pattern="*rep0.*\\.log$", full.names = TRUE)

for (i in seq(1,length(trees))){
  in_command <- " -b 20 -resample 100000 -log"
  for (j in seq(0,2)){
    in_command = paste(in_command, " ", gsub("rep0", paste("rep", j,sep=""), trees[i]), sep="")
  }
  
  out_command = gsub("rep0_", "", trees[i])
  out_command = gsub("individual", "combined", out_command)
  
  combined_command = gsub(".trees",".trees", out_command)
  combined_command = paste(" -o ", gsub("_rep0", "",combined_command), sep="")
  # combine the trees
  system(paste("/Applications/BEAST\\ 2.5.2/bin/logcombiner", in_command, gsub(".network.trees",".log", combined_command), sep=" "))
}
