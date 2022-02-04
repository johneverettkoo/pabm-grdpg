# simulation study for clustering algorithms
# compares OSC, SSC-ASE, SSC-A, and Louvain (stand-in for MM)
# simulations are for imbalanced community sizes
# configuration as-is parallelizes simulations across 10 threads
# which uses ~32 GB memory at its peak
# modify line 29 to increase/decrease parallelization (and memory use)

# packages, functions, etc.
import::from(magrittr, `%>%`)
import::from(foreach, foreach, `%do%`, `%dopar%`)
library(mclust)
library(ggplot2)
setwd('~/dev/pabm-grdpg')
source('simulations-and-data-analysis/functions.R')

# simulation parameters
K.vec <- c(2, 3, 4)
a1 <- 2
b1 <- 1
a2 <- 1
b2 <- 2
sparsity <- 1e-2
n.vec <- c(128, 256, 512, 1024, 2048, 4096)
iter <- 50
n.vec <- rev(n.vec)
K.vec <- rev(K.vec)
set.seed(314159)

doMC::registerDoMC(16)

# clustering simulation
clustering.df <- foreach(K = K.vec, .combine = dplyr::bind_rows) %do% {
  reorder.mat <- gtools::permutations(K, K)
  out.df <- foreach(n = n.vec, .combine = dplyr::bind_rows) %do% {
    print(paste('K =', K, ', n =', n))
    foreach(i = seq(iter), .combine = dplyr::bind_rows, .errorhandling = 'remove') %dopar% {
      Pz <- generate.P.beta(n, K, a1, a2, b1, b2, unbalanced = TRUE)
      P <- Pz$P
      z <- Pz$clustering
      A <- draw.graph(P)
      clustering.osc <- cluster.pabm(A, K, use.all = TRUE, normalize = FALSE)
      error.osc <- 1 - cluster.acc(clustering.osc, z, reorder.mat)
      clustering.ssc <- ssc(A, K, sparsity, normalize = TRUE, scale = FALSE,
                            offset = 1, d = K - 1)
      error.ssc <- 1 - cluster.acc(clustering.ssc, z, reorder.mat)
      clustering.louvain <- mod.max(A)
      error.louvain <- 1 - cluster.acc(clustering.louvain, z)
      clustering.ssc.A <- ssc2(A, K, sparsity, 
                               normalize = TRUE, parallel = FALSE)
      error.ssc.A <- 1 - cluster.acc(clustering.ssc.A, z)
      out <- dplyr::tibble(K = K, n = n,
                           error.osc = error.osc,
                           error.ssc.ase = error.ssc,
                           error.louvain = error.louvain,
                           error.ssc.A = error.ssc.A)
      # gc()
      return(out)
    } %>% 
      return()
  }
  gc()
  return(out.df)
}

gc()

# readr::write_csv(clustering.df, 'clustering-k-imbalanced.csv')

clustering.df %>%
  dplyr::group_by(n, K) %>%
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
  scale_x_log10(breaks = c(128, 256, 512, 1024, 2048, 4096),
                labels = c(expression(2^7), 
                           expression(2^8), 
                           expression(2^9), 
                           expression(2^10), 
                           expression(2^11), 
                           expression(2^12))) +
  theme(text = element_text(size = 20),
        legend.position = 'bottom') + 
  scale_y_log10() + 
  labs(y = 'error count', 
       colour = NULL, shape = NULL) +
  geom_line(aes(x = n, y = med.err * n,
                colour = 'OSC')) +
  geom_point(aes(x = n, y = med.err * n,
                 colour = 'OSC', shape = 'OSC'), size = 3) +
  geom_errorbar(aes(x = n, ymin = first.q * n, ymax = third.q * n,
                    colour = 'OSC'), width = .1) + 
  geom_line(aes(x = n, y = med.err.ssc * n,
                colour = 'SSC-ASE')) +
  geom_point(aes(x = n, y = med.err.ssc * n,
                 colour = 'SSC-ASE', shape = 'SSC-ASE'), size = 3) +
  geom_errorbar(aes(x = n, ymin = first.q.ssc * n, ymax = third.q.ssc * n,
                    colour = 'SSC-ASE'), width = .1) +
  geom_line(aes(x = n, y = med.err.ssc.A * n,
                colour = 'SSC-A')) +
  geom_point(aes(x = n, y = med.err.ssc.A * n,
                 colour = 'SSC-A', shape = 'SSC-A'), size = 3) +
  geom_errorbar(aes(x = n, ymin = first.q.ssc.A * n, ymax = third.q.ssc.A * n,
                    colour = 'SSC-A'), width = .1) +
  geom_line(aes(x = n, y = med.err.louvain * n,
                colour = 'MM-Louvain')) +
  geom_point(aes(x = n, y = med.err.louvain * n,
                 colour = 'MM-Louvain', shape = 'MM-Louvain'), size = 3) +
  geom_errorbar(aes(x = n, ymin = first.q.louvain * n, ymax = third.q.louvain * n,
                    colour = 'MM-Louvain'), width = .1) +
  scale_colour_brewer(palette = 'Set1') + 
  facet_wrap(~ K, labeller = 'label_both')
