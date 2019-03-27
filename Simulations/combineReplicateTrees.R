######################################################
######################################################
# combine the gene tree runs and run the mcc trees
######################################################
######################################################
library(ape)
# clear workspace
rm(list = ls())

# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)



trees <- list.files(path="./trees", pattern="rep1.*\\.trees$", full.names = TRUE)

system("rm -r combined")
system("mkdir combined")



# run log combiner
for (i in seq(1,length(trees))){
  in_command <- " -b 10 -resample 1000000 -log"
  for (j in seq(1,3)){
    in_command = paste(in_command, " ", gsub("rep1", paste("rep", j,sep=""), trees[i]), sep="")
  }
  
  out_command = gsub("_rep1", "", trees[i])
  out_command = gsub("/trees", "/combined", out_command)
  
  combined_command = gsub(".trees",".combined.trees", out_command)
  combined_command = paste(" -o ", combined_command, sep="")
  
  # combine the trees
  system(paste("/Applications/BEAST\\ 2.5.0/bin/logcombiner", in_command, combined_command, "", sep=" "))
  
}

