# packages, functions, etc.
import::from(magrittr, `%>%`)
import::from(foreach, foreach, `%do%`, `%dopar%`)
library(mclust)
library(ggplot2)
source('http://pages.iu.edu/~mtrosset/Courses/675/stress.r')
setwd('~/dev/pabm-grdpg')
import::here(embedding, draw.graph, draw.pabm.beta.2,
             normalized.laplacian,
             cluster.pabm, ssc, 
             embedding, 
             .from = 'functions.R')
library(ggplot2)

# parallel backend
doMC::registerDoMC(parallel::detectCores())

# simulation parameters
K <- 2
p <- K * (K + 1) / 2
q <- K * (K - 1) / 2
alpha <- .5
a1 <- 2
b1 <- 1
a2 <- 1
b2 <- 2
sparsity <- .01
# n.vec <- c(64, 128, 256, 512, 1024, 2048)
n.vec <- c(512, 1024, 2048, 3072, 4096)
iter <- 100
subsample <- 1e3
set.seed(314159)

# clustering simulation
densities.df <- foreach(n = n.vec, .combine = dplyr::bind_rows) %do% {
  Az <- draw.pabm.beta.2(n, a1, b1, a2, b2, alpha)
  A <- Az$A
  z <- Az$z
  A.eigen <- eigen(A, symmetric = TRUE)
  V <- A.eigen$vectors[, c(1:3, n)]
  B <- V %*% t(V)
  between.clust <- B[z == 1, z == 2]
  between.clust <- between.clust[upper.tri(between.clust)]
  between.clust <- sample(between.clust, subsample)
  dplyr::tibble(n = n, inner.prods = between.clust) %>% 
    return()
}

ggplot(densities.df) + 
  geom_density(aes(x = inner.prods * n, colour = factor(n)))

readr::write_csv(densities.df, 'densities.csv')
