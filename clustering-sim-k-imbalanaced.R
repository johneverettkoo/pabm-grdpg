# simulation using proposed clustering method
# run on IU RED (not a good idea to attempt on laptop)

# packages, functions, etc.
import::from(magrittr, `%>%`)
import::from(foreach, foreach, `%do%`, `%dopar%`)
library(mclust)
library(igraph)
library(ggplot2)
source('http://pages.iu.edu/~mtrosset/Courses/675/stress.r')
setwd('~/dev/pabm-grdpg')
import::here(embedding, draw.graph, 
             generate.P.beta,
             normalized.laplacian,
             cluster.pabm, ssc,
             .from = 'functions.R')

mod.max <- function(A) {
  Ag <- graph_from_adjacency_matrix(A, mode = 'undirected')
  clustering <- cluster_louvain(Ag)
  return(membership(clustering))
}

cluster.acc <- function(yhat, yobs) {
  table(yhat, yobs) %>% 
    as.matrix() %>% 
    apply(1, max) %>% 
    sum() %>% 
    magrittr::divide_by(length(yhat)) %>% 
    return()
}

# parallel backend
doMC::registerDoMC(18)

# simulation parameters
K.vec <- c(2, 3, 4)
a1 <- 2
b1 <- 1
a2 <- 1
b2 <- 2
sparsity <- 1e-3
n.vec <- c(128, 256, 512, 1024, 2048, 4096)
iter <- 50
# n.vec <- c(256, 512, 1024)
# iter <- 10
n.vec <- rev(n.vec)
K.vec <- rev(K.vec)
set.seed(314159)

# clustering simulation
clustering.df <- foreach(K = K.vec, .combine = dplyr::bind_rows) %do% {
  out.df <- foreach(n = n.vec, .combine = dplyr::bind_rows) %do% {
    print(paste('K =', K, ', n =', n))
    foreach(i = seq(iter), 
            .combine = dplyr::bind_rows, 
            .errorhandling = 'remove') %dopar% {
      Pz <- generate.P.beta(n, K, a1, a2, b1, b2, unbalanced = TRUE)
      P <- Pz$P
      z <- Pz$clustering
      A <- draw.graph(P)
      clustering <- cluster.pabm(A, K, use.all = TRUE)
      error <- 1 - cluster.acc(clustering, z)
      clustering.ssc <- ssc(A, K, sparsity)
      error.ssc <- 1 - cluster.acc(clustering.ssc, z)
      clustering.mm <- mod.max(A)
      error.mm <- 1 - cluster.acc(clustering.mm, z)
      dplyr::tibble(K = K, n = n, 
                    error = error, 
                    error.ssc = error.ssc,
                    error.mm = error.mm) %>% 
        return()
    } %>% 
      return()
  }
  gc()
  return(out.df)
}

gc()

clustering.df %>%
  dplyr::group_by(K, n) %>%
  dplyr::summarise(med.err = median(error),
                   first.q = quantile(error, .25),
                   third.q = quantile(error, .75),
                   med.err.ssc = median(error.ssc),
                   first.q.ssc = quantile(error.ssc, .25),
                   third.q.ssc = quantile(error.ssc, .75),
                   med.err.mm = median(error.mm),
                   first.q.mm = quantile(error.mm, .25),
                   third.q.mm = quantile(error.mm, .75)) %>%
  ggplot() +
  scale_y_log10() +
  # scale_x_log10() +
  labs(y = 'error') +
  geom_line(aes(x = n, y = med.err)) +
  geom_errorbar(aes(x = n, ymin = first.q, ymax = third.q)) +
  geom_line(aes(x = n, y = med.err.ssc, colour = 'ssc')) +
  geom_errorbar(aes(x = n, ymin = first.q.ssc, ymax = third.q.ssc,
                    colour = 'ssc')) + 
  geom_line(aes(x = n, y = med.err.mm, colour = 'mm')) + 
  geom_errorbar(aes(x = n, ymin = first.q.mm, ymax = third.q.mm,
                    colour = 'mm')) + 
  facet_wrap(~ K)

# export as csv
readr::write_csv(clustering.df, 'clustering-k-imbalanced.csv')
