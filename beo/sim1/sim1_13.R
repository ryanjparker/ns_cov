print(.libPaths())
.libPaths("/mnt/home/rjparker/Rlib")
print(.libPaths())
setwd("/mnt/home/rjparker/git/ns_cov")

which_exp <- 1
which_part <- 13
source("R/sim_study.R")
