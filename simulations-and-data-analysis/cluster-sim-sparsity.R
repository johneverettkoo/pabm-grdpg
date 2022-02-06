# simulation study for clustering algorithms
# compares OSC, SSC-ASE, SSC-A, and Louvain (stand-in for MM)
# for varying values of rho, the sparsity parameter
# simulations are for balanced community sizes
# configuration as-is parallelizes simulations across 16 threads

# packages, functions, etc.
import::from(magrittr, `%>%`)
import::from(foreach, foreach, `%do%`, `%dopar%`)
library(mclust)
library(ggplot2)
setwd('~/dev/pabm-grdpg')
source('simulations-and-data-analysis/functions.R')

doMC::registerDoMC(16)

# simulation parameters
n <- 2048
K <- 3
a1 <- 2
b1 <- 1
a2 <- 1
b2 <- 2
sparsity <- 1e-2  # for lasso within ssc
iter <- 50
rho.vec <- seq(.1, .9, .2)

# clustering simulation
clustering.df <- foreach(rho = rho.vec, .combine = dplyr::bind_rows) %do% {
  reorder.mat <- gtools::permutations(K, K)
  print(paste('rho =', rho))
  foreach(i = seq(iter), .combine = dplyr::bind_rows, .errorhandling = 'remove') %dopar% {
    Pz <- generate.P.beta(n, K, a1, a2, b1, b2)
    P <- Pz$P * rho
    z <- Pz$clustering
    A <- draw.graph(P)
    clustering.osc <- cluster.pabm(A, K, use.all = TRUE, normalize = FALSE)
    error.osc <- 1 - cluster.acc(clustering.osc, z, reorder.mat)
    clustering.ssc <- ssc(A, K, sparsity, normalize = TRUE, scale = FALSE)
    error.ssc <- 1 - cluster.acc(clustering.ssc, z, reorder.mat)
    clustering.louvain <- mod.max(A)
    error.louvain <- 1 - cluster.acc(clustering.louvain, z)
    clustering.ssc.A <- ssc2(A, K, sparsity, 
                             normalize = TRUE, parallel = FALSE)
    error.ssc.A <- 1 - cluster.acc(clustering.ssc.A, z)
    out.df <- dplyr::tibble(K = K, n = n, rho = rho, 
                            error.osc = error.osc,
                            error.ssc.ase = error.ssc,
                            error.louvain = error.louvain,
                            error.ssc.A = error.ssc.A)
    gc()
    return(out.df)
  } %>% 
    return()
}

gc()

readr::write_csv(clustering.df, 'clustering-sparsity.csv')

clustering.df %>%
  dplyr::group_by(rho, K, n) %>%
  dplyr::summarise(
    med.err = median(error.osc),
    first.q = quantile(error.osc, .25),
    third.q = quantile(error.osc, .75),
    med.err.ssc = median(error.ssc.ase),
    first.q.ssc = quantile(error.ssc.ase, .25),
    third.q.ssc = quantile(error.ssc.ase, .75),
    med.err.louvain = median(error.louvain),
    first.q.louvain = quantile(error.louvain, .25),
    third.q.louvain = quantile(error.louvain, .75),
    med.err.ssc.A = median(error.ssc.A),
    first.q.ssc.A = quantile(error.ssc.A, .25),
    third.q.ssc.A = quantile(error.ssc.A, .75)
  ) %>% 
  dplyr::ungroup() %>% 
  ggplot() +
  theme_bw() + 
  theme(text = element_text(size = 10)) + 
  labs(y = 'community detection error count', 
       colour = NULL, shape = NULL) +
  geom_line(aes(x = rho, y = med.err * n,
                colour = 'OSC')) +
  geom_point(aes(x = rho, y = med.err * n,
                 colour = 'OSC', shape = 'OSC')) +
  geom_errorbar(aes(x = rho, ymin = first.q * n, ymax = third.q * n,
                    colour = 'OSC'), width = .1) + 
  geom_line(aes(x = rho, y = med.err.ssc * n,
                colour = 'SSC-ASE')) +
  geom_point(aes(x = rho, y = med.err.ssc * n,
                 colour = 'SSC-ASE', shape = 'SSC-ASE')) +
  geom_errorbar(aes(x = rho, ymin = first.q.ssc * n, ymax = third.q.ssc * n,
                    colour = 'SSC-ASE'), width = .1) +
  geom_line(aes(x = rho, y = med.err.ssc.A * n,
                colour = 'SSC-A')) +
  geom_point(aes(x = rho, y = med.err.ssc.A * n,
                 colour = 'SSC-A', shape = 'SSC-A')) +
  geom_errorbar(aes(x = rho, ymin = first.q.ssc.A * n, ymax = third.q.ssc.A * n,
                    colour = 'SSC-A'), width = .1) +
  geom_line(aes(x = rho, y = med.err.louvain * n,
                colour = 'MM-Louvain')) + 
  geom_point(aes(x = rho, y = med.err.louvain * n,
                 colour = 'MM-Louvain', shape = 'MM-Louvain')) + 
  geom_errorbar(aes(x = rho, ymin = first.q.louvain * n, ymax = third.q.louvain * n,
                    colour = 'MM-Louvain'), width = .1) + 
  scale_colour_brewer(palette = 'Set1')
