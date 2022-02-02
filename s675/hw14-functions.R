sbm.edge.prob.matrix <- function(n1 = 100, 
                                 n2 = 100,
                                 p11 = .25,
                                 p22 = .25,
                                 p12 = .1) {
  # constructs an edge probability matrix for an SBM with K = 2
  # n1 and n2 are the number of members in communities 1 and 2
  # p11 is the within-community edge probability for community 1
  # p22 is the within-community edge probability for community 2
  # p12 is the between-community edge probability
  # returns a list with P (edge probability matrix) and z (community labels)
  
  P <- as.matrix(Matrix::bdiag(list(matrix(p11, n1, n1),
                                    matrix(p22, n2, n2))))
  P[P == 0] <- p12
  z <- c(rep(1, n1), rep(2, n2))
  return(list(P = P, z = z))
}

compute.error <- function(z, zhat) {
  # computes the error between true labels z and predicted labels zhat
  # this only works for K = 2
  
  error <- mean(z == zhat)
  if (error > .5) error <- 1 - error
  return(error)
}

draw.graph <- function(P) {
  # samples an adjacency matrix A from an edge probability matrix P
  
  n <- nrow(P)
  A <- matrix(0, nrow = n, ncol = n)
  A[upper.tri(A)] <- rbinom(n * (n - 1) / 2, 1, P[upper.tri(P)])
  A <- A + t(A)
  return(A)
}

ase <- function(A, p = 2, q = 0) {
  # adjacency spectral embedding
  
  n <- nrow(A)
  eigen.A <- eigen(A, symmetric = TRUE)
  
  if (p * q > 0) {
    keep <- c(seq(p), seq(n, n - q + 1))
  } else if (p == 0) {
    keep <- seq(n, n - q + 1)
  } else {
    keep <- seq(p)
  }
  
  U <- eigen.A$vectors[, keep]
  S <- diag(sqrt(abs(eigen.A$values[keep])))
  return(U %*% S)
}

estimate.sbm.probs.from.kmeans <- function(kmeans.out) {
  # estimates p11, p12, and p22 from kmeans clustering applied to the ASE of SBM
  
  centers <- kmeans.out$centers
  p11 <- as.numeric(crossprod(centers[1, ]))
  p22 <- as.numeric(crossprod(centers[2, ]))
  p12 <- as.numeric(crossprod(centers[1, ], centers[2, ]))
  
  return(list(p11 = p11, p12 = p12, p22 = p22))
}
