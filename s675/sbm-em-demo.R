library(mclust)
import::from('~/dev/pabm-grdpg/functions.R', draw.graph)

em.sbm <- function(A, eps = 1e-9, maxit = 100) {
  # EM for SBM with 2 communities
  # using mean field approximation
  
  # Pi is an n by 2 matrix of estimated community probabilities
  # Pi <- initial.guess(A)
  n <- nrow(A)
  Pi <- matrix(runif(n * 2), nrow = n, ncol = 2)
  Pi <- sweep(Pi, 1, rowSums(Pi), `/`)

  # Theta is a K by K matrix of community edge probabilities
  # start with some initial guess
  Theta <- matrix(c(.9, .1, .1, .9), 2, 2)
  
  d.Pi <- eps + 1
  niter <- 0
  
  while (d.Pi > eps) {
    Pi.old <- Pi
    
    Pi <- e.step(A, Pi, Theta)
    Theta <- m.step(A, Pi)
    
    niter <- niter + 1
    if (niter >= maxit) {
      warning('failed to converge')
      break
    }
    
    d.Pi <- norm(Pi - Pi.old, 'F') ** 2 / prod(dim(Pi))
    # print(d.Pi)
    zhat <- apply(Pi, 1, which.max)
    print(compute.error(z, zhat))
  }
  
  zhat <- apply(Pi, 1, which.max)
  
  return(list(Pi = Pi, Theta = Theta, z = zhat))
}

initial.guess <- function(A) {
  # GMM on the laplacian eigenmap
  
  n <- nrow(A)
  
  L <- diag(colSums(A)) - A
  L.eigen <- eigen(L, symmetric = TRUE)
  X <- L.eigen$vectors[, n - 1]
  Pi <- mclust::Mclust(X, 2)$z
  
  return(Pi)
}

e.step <- function(A, Pi, Theta) {
  n <- nrow(A)
  Pi.new <- matrix(0, nrow = n, ncol = 2)
  
  for (i in seq(n)) {
    for (k in seq(2)) {
      for (j in seq(n)[-i]) {
        for (l in seq(2)) {
          Pi.new[i, k] <- Pi.new[i, k] + 
            Pi[j, l] * 
            (A[i, j] * log(Theta[k, l]) + (1 - A[i, j]) * log(1 - Theta[k, l]))
        }
      }
    }
  }
  # Pi.new <- Pi.new - max(Pi.new)
  Pi.new <- exp(Pi.new)
  Pi.new <- sweep(Pi.new, 1, rowSums(Pi.new), `/`)
  return(Pi.new)
}

m.step <- function(A, Pi) {
  n <- nrow(A)
  Theta <- matrix(0, 2, 2)
  for (k in seq(2)) {
    for (l in seq(k, 2)) {
      numerator <- 0
      denominator <- 0
      for (j in seq(2, n)) {
        for (i in seq(1, j - 1)) {
          numerator <- numerator + A[i, j] * Pi[i, k] * Pi[j, l]
          denominator <- denominator + Pi[i, k] * Pi[j, l]
        }
      }
      Theta[k, l] <- numerator / denominator
    }
  }
  Theta[2, 1] <- Theta[1, 2]
  return(Theta)
}

compute.error <- function(z, zhat) {
  error <- mean(z == zhat)
  if (error > .5) error <- 1 - error
  return(error)
}

n1 <- 20
n2 <- 20
n <- n1 + n2
z <- c(rep(1, n1), rep(2, n2))
p.within <- .4
p.between <- .1
P <- matrix(p.between, nrow = n, ncol = n)
P[seq(n1), seq(n1)] <- p.within
P[seq(n1 + 1, n), seq(n1 + 1, n)] <- p.within
A <- draw.graph(P)
qgraph::qgraph(A, groups = factor(z))

Theta <- matrix(c(p.within, p.between, p.between, p.within), nrow = 2)
