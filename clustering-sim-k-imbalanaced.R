# simulation using proposed clustering method
# run on IU RED (not a good idea to attempt on laptop)

# packages, functions, etc.
import::from(magrittr, `%>%`)
import::from(foreach, foreach, `%do%`, `%dopar%`)
library(mclust)
library(igraph)
library(ggplot2)
source('http://pages.iu.edu/~mtrosset/Courses/675/stress.r')
setwd('~/dev/pabm-grdpg/codes_and_data')
source('EPalgosim.R')
source('functions.R')
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

cluster.sg <- function(A) {
  b.can = EPalgo(A,eps=0) # EP algorithm (no perturbation)
  Q.PA.can = rep(NA, ncol(b.can))	# array to store Q values
  Q.DC.can = rep(NA, ncol(b.can))	# array to store Q values
  for (i in 1:ncol(b.can)){
    #check if any cluster is empty
    foo = rep(NA, 2)
    for (clus in 1:2) {foo[clus]=sum(b.can[,i]==clus)}
    if (min(foo)==0) {stop('Empty groups are not allowed')} 
    Q.PA.can[i] = Q.PA(A, b=b.can[,i])   # fit PABM
    Q.DC.can[i] = Q.DC(A, b=b.can[,i])   # fit DCBM
  } # end of i for loop
  foo1 = order(-Q.PA.can)[1] 
  b.PA = b.can[,foo1]   # community assignment that maximises Q.PA
  return(b.PA)
}

# simulation parameters
K.vec <- c(2, 3, 4)
a1 <- 2
b1 <- 1
a2 <- 1
b2 <- 2
sparsity <- 1e-2
n.vec <- c(128, 256, 512, 1024, 2048, 4096)
iter <- 50
# n.vec <- c(256, 512, 1024)
# iter <- 10
n.vec <- rev(n.vec)
K.vec <- rev(K.vec)
cores.per.n <- 2 ** 12 * 20
set.seed(314159)

# clustering simulation
clustering.df <- foreach(K = K.vec, .combine = dplyr::bind_rows) %do% {
  out.df <- foreach(n = n.vec, .combine = dplyr::bind_rows) %do% {
    doMC::registerDoMC(min(parallel::detectCores(), 
                           floor(cores.per.n / n)))
    print(paste('K =', K, ', n =', n))
    foreach(i = seq(iter), 
            .combine = dplyr::bind_rows, 
            .errorhandling = 'remove') %dopar% {
      Pz <- generate.P.beta(n, K, a1, a2, b1, b2, unbalanced = TRUE)
      P <- Pz$P
      z <- Pz$clustering
      A <- draw.graph(P)
      clustering <- cluster.pabm(A, K, use.all = TRUE, normalize = FALSE)
      error <- 1 - cluster.acc(clustering, z)
      clustering.ssc <- ssc(A, K, K / n / 4, normalize = TRUE, scale = FALSE)
      error.ssc <- 1 - cluster.acc(clustering.ssc, z)
      if (K == 2) {
        clustering.mm <- cluster.sg(A)
      } else {
        clustering.mm <- mod.max(A)
      }
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
