print(.libPaths())
.libPaths("/mnt/home/rjparker/Rlib")
print(.libPaths())
setwd("/mnt/home/rjparker/git/ns_cov")

which_type <- WT
which_lambda <- WL
which_Nr <- WR
source("R/ozone.R")
