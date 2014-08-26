print(.libPaths())
.libPaths("/mnt/home/rjparker/Rlib")
print(.libPaths())
setwd("/mnt/home/rjparker/git/ns_cov")

which_exp <- 2
which_part <- 18
source("R/sim_study.R")
