# simulation using proposed clustering method
# run on IU RED (not a good idea to attempt on laptop)

# packages, functions, etc.
import::from(magrittr, `%>%`)
import::from(foreach, foreach, `%do%`, `%dopar%`)
library(mclust)
source('http://pages.iu.edu/~mtrosset/Courses/675/stress.r')
setwd('~/dev/pabm-grdpg')
import::here(embedding, draw.graph, draw.pabm.beta.2,
             normalized.laplacian,
             cluster.pabm, ssc,
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
sparsity <- .01
n.vec <- c(64, 128, 256, 512, 1024, 2048)
iter <- 100
set.seed(314159)

# clustering simulation
clustering.df <- foreach(n = n.vec, .combine = dplyr::bind_rows) %do% {
  foreach(i = seq(iter), 
          .combine = dplyr::bind_rows, 
          .errorhandling = 'remove') %dopar% {
    Az <- draw.pabm.beta.2(n, a1, b1, a2, b2, alpha)
    A <- Az$A
    z <- Az$z
    clustering <- cluster.pabm(A, K)
    error <- mean(clustering == z)
    if (error > .5) error <- 1 - error
    clustering.ssc <- ssc(A, K, sparsity)
    error.ssc <- mean(clustering.ssc == z)
    if (error.ssc > .5) error.ssc <- 1 - error.ssc
    dplyr::tibble(n = n, error = error, error.ssc = error.ssc) %>% 
      return()
  } %>% 
    return()
}

# clustering.df %>%
#   dplyr::group_by(n) %>%
#   dplyr::summarise(med.err = median(error),
#                    first.q = quantile(error, .25),
#                    third.q = quantile(error, .75),
#                    med.err.ssc = median(error.ssc),
#                    first.q.ssc = quantile(error.ssc, .25),
#                    third.q.ssc = quantile(error.ssc, .75)) %>%
#   ggplot() +
#   # scale_y_log10() +
#   # scale_x_log10() +
#   labs(y = 'error') +
#   geom_line(aes(x = n, y = med.err)) +
#   geom_errorbar(aes(x = n, ymin = first.q, ymax = third.q)) +
#   geom_line(aes(x = n, y = med.err.ssc, colour = 'ssc')) +
#   geom_errorbar(aes(x = n, ymin = first.q.ssc, ymax = third.q.ssc,
#                     colour = 'ssc'))

# export as csv
readr::write_csv(clustering.df, 'clustering.csv')
