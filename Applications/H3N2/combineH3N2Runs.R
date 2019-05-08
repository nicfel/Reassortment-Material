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


system("rm -r combined")
system("mkdir combined")


trees <- list.files(path=paste("./out/",sep=""), pattern="*rep0.*rk\\.trees$", full.names = TRUE)

for (i in seq(1,length(trees))){
  in_command <- " -b 20 -log"
  for (j in seq(0,2)){
    in_command = paste(in_command, " ", gsub("rep0", paste("rep", j,sep=""), trees[i]), sep="")
  }
  
  out_command = gsub("rep0_", "", trees[i])
  out_command = gsub("out", "combined", out_command)
  
  combined_command = gsub(".trees",".trees", out_command)
  combined_command = paste(" -o ", gsub("_rep0", "",combined_command), sep="")
  # combine the trees
  system(paste("/Applications/BEAST\\ 2.5.2/bin/logcombiner", in_command, combined_command, "", sep=" "))
  system(paste("/Applications/BEAST\\ 2.5.2/bin/logcombiner", gsub(".network.trees",".log", in_command), gsub(".network.trees",".log", combined_command), sep=" "))
}

# compute how many reassortment events happen on the trunk vs. not on the trunk
networks <- list.files(path=paste("./combined/",sep=""), pattern="*rk\\.trees$", full.names = TRUE)
for (i in seq(1,length(networks))){
  system(paste("java -jar ./../../Software/TrunkReassortment.jar -burnin 0 -removeSegments 8 -trunkDefinition minTipDistance -minTipDistance 4",
               networks[[i]], gsub("network.trees", "trunk.txt", networks[[i]])))
  system(paste("java -jar ./../../Software/ReassortmentDistance.jar -burnin 90",
               networks[[i]], gsub("network.trees", "distance.txt", networks[[i]])))
  system(paste("java -jar ./../../Software/ReassortmentNetworkSummarizer.jar -burnin 0 -removeSegments 8",
               networks[[i]], gsub("network.trees", "summary.trees", networks[[i]])))
  
}



# for the dummy networks, also calculate the pairwise reassortment distances

# system("java -jar ./../../Software/TrunkReassortment.jar ./combined/h3n2_rep0.network.trees combined/trunk_reassortment.txt")
# system("java -jar ./../../Software/TrunkReassortment.jar ./combined/h3n2old_rep0.network.trees combined/trunk_old_reassortment.txt")
# system("java -jar ./../../Software/TrunkReassortment.jar ./combined/h3n2recent_rep0.network.trees combined/trunk_recent_reassortment.txt")


# # simulate reassortment networks from the prior with estimated ne and reassortment rates
# system("java -jar ./../../Software/Coevo.jar ./simulation/h3n2_sim.xml")
# system("java -jar ./../../Software/Coevo.jar ./simulation/h3n2old_sim.xml")
# system("java -jar ./../../Software/Coevo.jar ./simulation/h3n2recent_sim.xml")
# 
# 
# system("java -jar ./../../Software/TrunkReassortment.jar ./simulation/h3n2_sim.network.trees simulation/trunk_reassortment.txt")
# system("java -jar ./../../Software/TrunkReassortment.jar ./simulation/h3n2old_sim.network.trees simulation/trunk_old_reassortment.txt")
# system("java -jar ./../../Software/TrunkReassortment.jar ./simulation/h3n2recent_sim.network.trees simulation/trunk_recent_reassortment.txt")
