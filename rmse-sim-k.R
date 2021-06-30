# simulation using proposed parameter estimation method
# run on IU RED (not a good idea to attempt on laptop)

# packgaes, functions, etc.
import::from(magrittr, `%>%`)
import::from(foreach, foreach, `%do%`)
library(ggplot2)
setwd('~/dev/pabm-grdpg')
import::here(embedding, draw.graph, I.pq, generate.P.beta,
             estimate.lambda.block, lambda.rmse, lambda.rmse.mle,
             .from = 'functions.R')

# parallel backend
doMC::registerDoMC(parallel::detectCores() / 2)

# simulation params
K.vec <- c(4, 3, 2)
a1 <- b2 <- 2
a2 <- b1 <- 1
n.vec <- c(128, 256, 512, 1024, 2048, 4096)
n.vec <- rev(n.vec)
iter <- 50
set.seed(314159)

# simulation
rmse.df <- foreach(K = K.vec, .combine = dplyr::bind_rows) %do% {
  p <- K * (K + 1) / 2
  q <- K * (K - 1) / 2
  Ipq <- I.pq(p, q)
  
  out <- foreach(n = n.vec, .combine = dplyr::bind_rows) %do% {
    print(paste('n =', n, 'K =', K))
    plyr::ldply(seq(iter), function(i) {
      P.z <- generate.P.beta(n, K)
      P <- P.z$P
      z <- P.z$clustering
      A <- draw.graph(P)
      rmse <- lambda.rmse(P, A, z)
      rmse.mle <- lambda.rmse.mle(P, A, z)
      dplyr::tibble(K = K, n = n, rmse = rmse, rmse.mle = rmse.mle) %>% 
        return()
    }, .parallel = TRUE) %>% 
      return()
  }
  gc()
  return(out)
}

rmse.df %>%
  na.omit() %>%
  dplyr::group_by(K, n) %>%
  dplyr::summarise(median.rmse = median(rmse),
                   q1.rmse = quantile(rmse, .25),
                   q3.rmse = quantile(rmse, .75),
                   median.rmse.mle = median(rmse.mle),
                   q1.rmse.mle = quantile(rmse.mle, .25),
                   q3.rmse.mle = quantile(rmse.mle, .75)) %>%
  dplyr::ungroup() %>%
  ggplot() +
  scale_y_log10() +
  scale_x_log10() +
  geom_line(aes(x = n, y = median.rmse, colour = 'Proposed')) +
  geom_errorbar(aes(x = n, ymin = q1.rmse, ymax = q3.rmse,
                    colour = 'Proposed')) +
  geom_line(aes(x = n, y = median.rmse.mle, colour = 'MLE-based')) +
  geom_errorbar(aes(x = n, ymin = q1.rmse.mle, ymax = q3.rmse.mle,
                    colour = 'MLE-based')) +
  scale_colour_brewer(palette = 'Set1') +
  labs(y = 'RMSE', colour = NULL) +
  facet_wrap(~ K)

readr::write_csv(rmse.df, 'rmse-k.csv')
