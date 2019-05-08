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


system("rm -r trees/combined")
system("mkdir trees/combined")


system("rm -r trees/mcc")
system("mkdir trees/mcc")


# define the segments
segments = c("HA", "MP", "NA", "NP", "NS1", "PA", "PB1", "PB2")


trees <- list.files(path=paste("./trees/individual",sep=""), pattern="*rep0.*\\.trees$", full.names = TRUE)

for (i in seq(1,length(trees))){
  in_command <- " -b 40 -log"
  for (j in seq(0,2)){
    in_command = paste(in_command, " ", gsub("rep0", paste("rep", j,sep=""), trees[i]), sep="")
  }

  out_command = gsub("rep0_", "", trees[i])
  out_command = gsub("individual", "combined", out_command)

  combined_command = paste(" -o ", gsub("_rep0", "",out_command), sep="")
  # combine the trees
  system(paste("/Applications/BEAST\\ 2.5.2/bin/logcombiner", in_command, gsub("_rep0", "", combined_command), "", sep=" "))
  tree_name = gsub("-o ", "", combined_command)
  system(paste("/Applications/BEAST\\ 2.5.2/bin/treeannotator -burnin 0", tree_name, gsub("/combined/", "/mcc/", tree_name)))

}


# compute the clade credibilities between accounting for reassortment and not accounting for it
system("rm -r trees/clades")
system("mkdir trees/clades")

timewindow = c("" ,"new")
for (i in seq(1,length(timewindow))){
  for (j in seq(1, length(segments))){
    target = paste("./trees/mcc/h3n2", timewindow[[i]], ".", segments[[j]], ".trees", sep="")
    tree = paste("./trees/combined/h3n2", timewindow[[i]], "norea.", segments[[j]], ".trees", sep="")
    print(tree)
    outfile = paste("./trees/clades/h3n2", timewindow[[i]], ".", segments[[j]], ".trees", sep="")
  
    system(paste("/Applications/BEAST\\ 2.5.2/bin/treeannotator -heights mean -burnin 0 -target",target, tree,  outfile))
  
  }
}
