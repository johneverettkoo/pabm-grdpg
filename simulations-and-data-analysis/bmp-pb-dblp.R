# british mp, political blogs, and dblp data analysis
# original data and code for data processing found at 
# https://sites.google.com/vt.edu/sengupta/software-and-datasets
# included data have already been transformed using the above scripts
# not including analysis using modularity maximization
# see paper by sengupta and chen for their results

library(mclust)
import::from(magrittr, `%>%`)
setwd('~/dev/pabm-grdpg')
source('simulations-and-data-analysis/functions.R')

doMC::registerDoMC(8)

compute.errors <- function(A, z, K = 2, 
                           normalize = TRUE, use.all = TRUE,
                           p = K * (K + 1) / 2, q = K * (K - 1) / 2, d.osc = K + 1, 
                           lambda = 1e-1, scale = TRUE, normalize.ssc = TRUE,
                           d = K, offset = 0,
                           parallel = TRUE) {
  zhat.osc <- cluster.pabm(A, K, 
                           normalize = normalize, use.all = use.all,
                           p = p, q = q, d.eigenmap = d.osc)
  zhat.ssc <- ssc(A, K, lambda, scale = scale, normalize = normalize.ssc,
                  d = d, offset = offset,
                  parallel = parallel)
  zhat.mm <- cluster.sg(A)
  list(error.osc = 1 - cluster.acc(z, zhat.osc),
       error.ssc = 1 - cluster.acc(z, zhat.ssc),
       error.mm = 1 - cluster.acc(z, zhat.mm))
}

bmp <- readRDS('simulations-and-data-analysis/bmp.rds')
pb <- readRDS('simulations-and-data-analysis/pb.rds')
dblp <- readRDS('simulations-and-data-analysis/dblp.rds')

compute.errors(bmp$A, bmp$z, 2)
compute.errors(pb$A, pb$z, 2, 
               normalize = FALSE, p = 2, q = 0, d.osc = 2, 
               lambda = .005, normalize.ssc = FALSE, scale = TRUE, 
               offset = 0, d = 2)
compute.errors(dblp$A, dblp$z, 2, 
               normalize = FALSE, 
               lambda = 1e-1, normalize.ssc = FALSE, scale = FALSE,
               d = 2, offset = 1)
