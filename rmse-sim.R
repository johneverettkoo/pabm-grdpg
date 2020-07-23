# simulation using proposed parameter estimation method
# run on IU RED (not a good idea to attempt on laptop)

# packages, functions, etc.
import::from(magrittr, `%>%`)
import::from(foreach, foreach, `%do%`)
setwd('~/dev/pabm-grdpg')
import::here(embedding, draw.graph, I.pq, 
             estimate.lambda.block, lambda.rmse, lambda.rmse.mle,
             .from = 'functions.R')

# parallel backend
doMC::registerDoMC(parallel::detectCores())

# simulation params
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
Ipq <- I.pq(p, q)
set.seed(314159)

# simulation
rmse.df <- foreach(n = n.vec, .combine = dplyr::bind_rows) %do% {
  plyr::ldply(seq(iter), function(i) {
    z <- sample(c(1, 2), n, prob = c(alpha, 1 - alpha), replace = TRUE)
    z <- sort(z)
    
    n1 <- sum(z == 1)
    n2 <- sum(z == 2)
    
    lambda11 <- rbeta(n1, a1, b1)
    lambda22 <- rbeta(n2, a1, b1)
    lambda12 <- rbeta(n1, a2, b2)
    lambda21 <- rbeta(n2, a2, b2)
    
    X <- cbind(c(lambda11, rep(0, n2)),
               c(lambda12, rep(0, n2)),
               c(rep(0, n1), lambda21),
               c(rep(0, n1), lambda22))
    Y <- cbind(c(lambda11, rep(0, n2)),
               c(rep(0, n1), lambda21),
               c(lambda12, rep(0, n2)),
               c(rep(0, n1), lambda22))
    P <- X %*% t(Y)
    A <- draw.graph(P)
    lambda.matrix <- cbind(c(lambda11, lambda21),
                           c(lambda12, lambda22))
    Zhat <- embedding(A, p, q)
    Phat <- Zhat %*% Ipq %*% t(Zhat)
    rmse <- lambda.rmse(lambda.matrix, Phat, z)
    rmse.mle <- lambda.rmse.mle(lambda.matrix, A, z)
    dplyr::tibble(n = n, rmse = rmse, rmse.mle = rmse.mle) %>% 
      return()
  }, .parallel = TRUE) %>% 
    return()
}

readr::write_csv(rmse.df, 'rmse.csv')
