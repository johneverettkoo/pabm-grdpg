# simulation using proposed clustering method
# run on IU RED (not a good idea to attempt on laptop)

# packages, functions, etc.
import::from(magrittr, `%>%`)
import::from(foreach, foreach, `%do%`, `%dopar%`)
library(mclust)
source('http://pages.iu.edu/~mtrosset/Courses/675/stress.r')
setwd('~/dev/pabm-grdpg')
import::here(embedding, draw.graph, 
             generate.P.beta,
             normalized.laplacian,
             cluster.pabm, ssc,
             .from = 'functions.R')

# parallel backend
doMC::registerDoMC(parallel::detectCores())

# simulation parameters
K.vec <- c(2, 4, 8)
a1 <- 2
b1 <- 1
a2 <- 1
b2 <- 2
sparsity <- .01
n.vec <- c(64, 128, 256, 512, 1024, 2048)
iter <- 100
set.seed(314159)

# clustering simulation
clustering.df <- foreach(K = K.vec, .combine = dplyr::bind_rows) %do% {
  foreach(n = n.vec, .combine = dplyr::bind_rows) %do% {
    print(paste('K =', K, ', n =', n))
    foreach(i = seq(iter), 
            .combine = dplyr::bind_rows, 
            .errorhandling = 'remove') %dopar% {
      Pz <- generate.P.beta(n, K, a1, a2, b1, b2)
      P <- Pz$P
      z <- Pz$clustering
      A <- draw.graph(P)
      clustering <- cluster.pabm(A, K)
      error <- 1 - fossil::adj.rand.index(z, clustering)
      clustering.ssc <- ssc(A, K, sparsity)
      error.ssc <- 1 - fossil::adj.rand.index(z, clustering.ssc)
      dplyr::tibble(K = K, n = n, error = error, 
                    error.ssc = error.ssc) %>% 
        return()
    } %>% 
      return()
  } %>%
    return()
}

clustering.df %>%
  dplyr::group_by(K, n) %>%
  dplyr::summarise(med.err = median(error),
                   first.q = quantile(error, .25),
                   third.q = quantile(error, .75),
                   med.err.ssc = median(error.ssc),
                   first.q.ssc = quantile(error.ssc, .25),
                   third.q.ssc = quantile(error.ssc, .75)) %>%
  ggplot() +
  # scale_y_log10() +
  # scale_x_log10() +
  labs(y = 'error') +
  geom_line(aes(x = n, y = med.err)) +
  geom_errorbar(aes(x = n, ymin = first.q, ymax = third.q)) +
  geom_line(aes(x = n, y = med.err.ssc, colour = 'ssc')) +
  geom_errorbar(aes(x = n, ymin = first.q.ssc, ymax = third.q.ssc,
                    colour = 'ssc')) + 
  facet_wrap(~ K)

# export as csv
readr::write_csv(clustering.df, 'clustering-k.csv')
