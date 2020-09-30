# simulation using proposed clustering method
# run on IU RED (not a good idea to attempt on laptop)
# this code is mostly from scripts provided by sengupta & chen

# packages, functions, etc.
import::from(magrittr, `%>%`)
import::from(foreach, foreach, `%do%`, `%dopar%`)
library(mclust)
library(ggplot2)
library(igraph)
library(clue)
source('http://pages.iu.edu/~mtrosset/Courses/675/stress.r')
setwd('~/dev/pabm-grdpg')
import::here(embedding, draw.graph, 
             generate.P.beta,
             normalized.laplacian,
             cluster.pabm, ssc,
             .from = 'functions.R')
setwd('codes_and_data')
source('functions.R')
source('EPalgosim.R')
setwd('..')

mod.max.clust <- function(A) {
  # code copy-pasted from sengupta & chen
  Eps <- (0:19)/20
  b.can <- NULL
  for (ieps in 1:length(Eps)){
    b.can<-c(b.can,list(EPalgo(A,eps=Eps[ieps]))) # eps=0 for unperturbed version
  }
  
  b.can.unpert <- b.can[[1]]     # matrix of candidate assignments for eps = 0
  n.can.unpert <- dim(b.can.unpert)[2]  # total number of candidate assignments
  Q.PA.unpert = rep(NA, n.can.unpert)  # array to store modularity values
  
  ##### modularity calculation #####
  ptm<-proc.time()
  for (ican in 1:n.can.unpert){
    Q.PA.unpert[ican] = Q.PA(A, b=b.can.unpert[,ican])
  }
  index.PA<-which(Q.PA.unpert==max(Q.PA.unpert))
  b.PA.unpert = b.can.unpert[,index.PA]
  
  return(b.PA.unpert)
}

# parallel backend
doMC::registerDoMC(parallel::detectCores())

# simulation parameters
K <- 2
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
  print(paste('K =', K, ', n =', n))
  foreach(i = seq(iter), 
          .combine = dplyr::bind_rows, 
          .errorhandling = 'remove') %dopar% {
            Pz <- generate.P.beta(n, K, a1, a2, b1, b2)
            P <- Pz$P
            z <- Pz$clustering
            A <- draw.graph(P)
            clustering <- mod.max.clust(A)
            error <- 1 - fossil::adj.rand.index(z, clustering)
            # P.hat <- P_PA(A, z)
            # mu.PA.hat <- mu.hat(P.hat,z)
            dplyr::tibble(K = K, n = n, error.mm = error) %>% 
              return()
          } %>% 
    return()
}

clustering.df %>%
  dplyr::group_by(K, n) %>%
  dplyr::summarise(med.err = median(error.mm),
                   first.q = quantile(error.mm, .25),
                   third.q = quantile(error.mm, .75)) %>%
  ggplot() +
  scale_y_log10() +
  # scale_x_log10() +
  labs(y = 'error') +
  geom_line(aes(x = n, y = med.err)) +
  geom_errorbar(aes(x = n, ymin = first.q, ymax = third.q)) +
  facet_wrap(~ K)

# export as csv
readr::write_csv(clustering.df, 'mod-max-k.csv')
