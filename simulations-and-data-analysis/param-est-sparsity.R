# simulation using proposed parameter estimation method
# for varying rho, the sparsity parameter
# compared against mle-based method from sengupta and chen

# packgaes, functions, etc.
import::from(magrittr, `%>%`)
import::from(foreach, foreach, `%do%`)
library(ggplot2)
setwd('~/dev/pabm-grdpg')
source('simulations-and-data-analysis/functions.R')

doMC::registerDoMC(16)

# simulation params
K <- 3
a1 <- b2 <- 2
a2 <- b1 <- 1
n <- 2048
rho.vec <- seq(.1, 1, .1)
iter <- 50

# simulation
rmse.df <- foreach(rho = rho.vec, .combine = dplyr::bind_rows) %do% {
  print(paste('rho =', rho))
  plyr::ldply(seq(iter), function(i) {
    P.z <- generate.P.beta(n, K,
                           a1 = a1, a2 = a2, 
                           b1 = b1, b2 = b2)
    P <- P.z$P * rho
    z <- P.z$clustering
    A <- draw.graph(P)
    rmse <- lambda.rmse(P, A, z, rho)
    rmse.mle <- lambda.rmse.mle(P, A, z, rho)
    dplyr::tibble(K = K, n = n, rho = rho, 
                  rmse = rmse, rmse.mle = rmse.mle) %>%
      return()
  }, .parallel = TRUE) %>%
    return()
}

readr::write_csv(rmse.df, 'rmse-sparsity.csv')

rmse.df %>% 
  na.omit() %>% 
  dplyr::group_by(rho) %>% 
  dplyr::summarise(median.rmse = median(rmse),
                   q1.rmse = quantile(rmse, .25),
                   q3.rmse = quantile(rmse, .75),
                   median.rmse.mle = median(rmse.mle),
                   q1.rmse.mle = quantile(rmse.mle, .25),
                   q3.rmse.mle = quantile(rmse.mle, .75)) %>% 
  dplyr::ungroup() %>% 
  ggplot() +
  geom_line(aes(x = rho, y = median.rmse, colour = 'Algorithm 3')) + 
  geom_point(aes(x = rho, y = median.rmse, colour = 'Algorithm 3', 
                 shape = 'Algorithm 3'), size = 3) + 
  geom_errorbar(aes(x = rho, ymin = q1.rmse, ymax = q3.rmse,
                    colour = 'Algorithm 3'), width = .1) + 
  geom_line(aes(x = rho, y = median.rmse.mle, colour = 'MLE-based')) + 
  geom_point(aes(x = rho, y = median.rmse.mle, colour = 'MLE-based',
                 shape = 'MLE-based'), size = 3) + 
  geom_errorbar(aes(x = rho, ymin = q1.rmse.mle, ymax = q3.rmse.mle,
                    colour = 'MLE-based'), width = .1) + 
  scale_colour_brewer(palette = 'Set1') + 
  labs(x = expression(rho), y = 'RMSE', colour = NULL, shape = NULL) + 
  scale_y_log10() + 
  theme_bw()
