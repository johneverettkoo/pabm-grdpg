# simulation using proposed clustering method
# looking specifically at distributions of between-cluster inner products
# run on IU RED (not a good idea to attempt on laptop)

# packages, functions, etc.
import::from(magrittr, `%>%`)
import::from(foreach, foreach, `%do%`)
library(mclust)
source('http://pages.iu.edu/~mtrosset/Courses/675/stress.r')
setwd('~/dev/pabm-grdpg')
import::here(embedding, draw.graph, draw.pabm.beta.2,
             normalized.laplacian,
             cluster.pabm,
             .from = 'functions.R')
# library(ggplot2)

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
n.vec <- c(64, 128, 256, 512, 1024)
iter <- 100
set.seed(314159)

# simulation
densities.df <- plyr::ldply(n.vec, function(n) {
  Az <- draw.pabm.beta.2(n, a1, b1, a2, b2, alpha)
  A <- Az$A
  z <- Az$z
  n1 <- sum(z == 1)
  n2 <- sum(z == 2)
  V <- eigen(A, symmetric = TRUE)$vectors[, c(seq(3), n)]
  B <- V %*% t(V)
  inner.prods <- as.numeric(B[(n1+1):n, (n1+1):n])
  return(dplyr::tibble(n = n, inner.prods = inner.prods))
}, .parallel = TRUE)

readr::write_csv(densities.df, 'densities.csv')
