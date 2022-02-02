# leeds butterfly data
# the data is truncated and transformed based on descriptions from noroozi, rimal, pensky
# 

library(ggplot2)
library(mclust)
import::from(magrittr, `%>%`)
setwd('~/dev/pabm-grdpg')
source('simulations-and-data-analysis/functions.R')

n.edges <- 20566  # based on description from noroozi, rimal, pensky

# data from https://github.com/wangboyunze/Network_Enhancement
butterfly <- R.matlab::readMat('data/Raw_butterfly_network.mat')

A <- butterfly$W.butterfly0
z <- butterfly$labels

# based on description from noroozi, rimal, pensky
keep <- c(6, 2, 9, 4)
A <- A[z %in% keep, z %in% keep]
z <- z[z %in% keep]
n <- length(z)
K <- 4
delta <- quantile(A[upper.tri(A)], 1 - n.edges / (n * (n - 1) / 2))
A <- (A > delta) * 1

# OSC
z.hat <- cluster.pabm(A, K)
ari <- round(fossil::adj.rand.index(z, z.hat) * 100)
plot.A(A, z.hat)

# SSC-ASE
z.hat.ssc <- ssc(A, K, 1e-1)
ari.ssc <- round(fossil::adj.rand.index(z, z.hat.ssc) * 100)
