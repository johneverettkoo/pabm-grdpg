embedding <- function(A, p = NULL, q = NULL,
                      eps = 1e-6) {
  eigen.A <- eigen(A, symmetric = TRUE)
  if (is.null(p) | is.null(q)) {
    keep <- (abs(eigen.A$values) > eps)
  } else {
    keep <- c(seq(p), seq(n, n - q + 1))
  }
  
  U <- eigen.A$vectors[, keep]
  S <- diag(sqrt(abs(eigen.A$values[keep])))
  return(U %*% S)
}

draw.graph <- function(P) {
  A <- apply(P, 1:2, function(p) rbinom(1, 1, p))
  A[lower.tri(A)] <- 0
  diag(A) <- 0
  A <- A + t(A)
  return(A)
}

I.pq <- function(p, q) {
  return(diag(c(rep(1, p), rep(-1, q))))
}

draw.pabm.beta.2 <- function(n, a1, b1, a2, b2, alpha = .5) {
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
  
  return(list(A = A, z = z))
}

normalized.laplacian <- function(W) {
  n <- nrow(W)
  d.neg.sqrt <- colSums(W) ** -.5
  d.neg.sqrt[is.infinite(d.neg.sqrt)] <- 0
  D.neg.sqrt <- diag(d.neg.sqrt)
  I <- diag(n)
  return(I - D.neg.sqrt %*% W %*% D.neg.sqrt)
}

cluster.pabm <- function(A, K) {
  p <- K * (K + 1) / 2
  q <- K * (K - 1) / 2
  V <- eigen(A, symmetric = TRUE)$vectors[, c(seq(p), seq(n, n - q + 1))]
  V <- V / mean(V)
  B <- abs(V %*% t(V))
  # L <- graph.laplacian(B)
  L <- normalized.laplacian(B)
  eigenmap <- eigen(L, symmetric = TRUE)$vectors[, seq(n - 1, n - K)]
  clustering <- mclust::Mclust(eigenmap, K)$classification
}

rot.mat.2 <- function(angle, axis = 1) {
  # rotation matrix for K=2
  
  if (axis %in% c(1, 3)) {
    rotation <- matrix(c(cos(angle), -sin(angle), sin(angle), cos(angle)),
                       nrow = 2, ncol = 2, byrow = TRUE)
    if (axis == 1) {
      return(as.matrix(Matrix::bdiag(1, rotation, -1)))
    } else {
      return(as.matrix(Matrix::bdiag(rotation, 1, -1)))
    }
  } else if (axis == 2) {
    return(matrix(c(cos(angle), 0, sin(angle), 0,
                    0, 1, 0, 0,
                    -sin(angle), 0, cos(angle), 0,
                    0, 0, 0, -1),
                  nrow = 4, ncol = 4,
                  byrow = TRUE))
  } else {
    stop(simpleError('axis must be 1, 2, or 3'))
  }
}

hyp.rot.mat.2 <- function(angle, axis = 1) {
  # hyperbolic rotation matrix for $K=2
  
  if (axis == 1) {
    return(matrix(c(1, 0, 0, 0,
                    0, 1, 0, 0, 
                    0, 0, cosh(angle), sinh(angle),
                    0, 0, sinh(angle), cosh(angle)), 
                  byrow = TRUE,
                  nrow = 4, ncol = 4))
  } else if (axis == 2) {
    return(matrix(c(1, 0, 0, 0,
                    0, cosh(angle), 0, sinh(angle),
                    0, 0, 1, 0,
                    0, sinh(angle), 0, cosh(angle)),
                  byrow = TRUE,
                  nrow = 4, ncol = 4))
  } else if (axis == 3) {
    return(matrix(c(cosh(angle), 0, 0, sinh(angle),
                    0, 1, 0, 0, 
                    0, 0, 1, 0, 
                    sinh(angle), 0, 0, cosh(angle)),
                  byrow = TRUE,
                  nrow = 4, ncol = 4))
  } else {
    stop(simpleError('axis must be 1, 2, or 3'))
  }
}

construct.Q <- function(angles) {
  rot.angles <- angles[1:3]
  hyp.angles <- c(angles[1], angles[4], -angles[4])
  rot.matrices <- lapply(seq(3), function(i) rot.mat.2(rot.angles[i], i))
  hyp.matrices <- lapply(seq(3), function(i) hyp.rot.mat.2(hyp.angles[i], i))
  matrices <- c(rot.matrices, hyp.matrices)
  return(Reduce(`%*%`, matrices))
}

obj.fun <- function(angles, Z.hat, clusters) {
  R <- construct.Q(angles)
  
  XU.c <- Z.hat %*% R
  in.prods <- XU.c[clusters == 1, ] %*% t(XU.c[clusters == 2, , drop = FALSE])
  return(
    sum(XU.c[clusters == 1, 2] ** 2) + 
      sum(XU.c[clusters == 2, 1] ** 2) + 
      sum((XU.c[clusters == 1, 3] - XU.c[clusters == 1, 4]) ** 2) + 
      sum((XU.c[clusters == 2, 3] + XU.c[clusters == 2, 4]) ** 2)
  )
}

transform.embedding <- function(Z.hat, clusters, init.angles) {
  # init.angles <- runif(4, -pi, pi)
  # init.angles <- rep(0, 4)
  angles.1 <- optim(init.angles, function(x) obj.fun(x, Z.hat, clusters), 
                    method = 'SANN')$par
  angles.2 <- optim(angles.1, function(x) obj.fun(x, Z.hat, clusters),
                    method = 'Nelder-Mead',
                    control = list(abstol = 1e-10, maxit = 1e3))$par
  R <- construct.Q(angles.2)
  return(Z.hat %*% R)
}

estimate.lambda.block <- function(P.block, within = FALSE) {
  if (within) {
    P.svd <- svd(P.block)
    u <- P.svd$u[, 1]
    if (mean(u) < 0) u <- -u
    s <- sqrt(P.svd$d[1])
    return(s * u)
  } else {
    P.svd <- svd(P.block)
    u <- P.svd$u[, 1]
    if (mean(u) < 0) u <- -u
    v <- P.svd$v[, 1]
    if (mean(v) < 0) v <- -v
    s <- sqrt(P.svd$d[1])
    return(list(s * u, s * v))
  }
}

lambda.rmse <- function(lambda.matrix, Phat, clustering) {
  K <- max(clustering)
  n <- nrow(Phat)
  n.vector <- sapply(seq(K), function(k) sum(clustering == k))
  sapply(seq(K), function(k) {
    low.ind.k <- ifelse(k == 1, 1, sum(n.vector[seq(k - 1)]) + 1)
    high.ind.k <- sum(n.vector[seq(k)])
    sapply(seq(k), function(l) {
      P.block <- Phat[clustering == k, clustering == l]
      if (k == l) {
        lambda <- as.numeric(lambda.matrix[seq(low.ind.k, high.ind.k), l])
        lambda.hat <- estimate.lambda.block(P.block, TRUE)
        return(sum((lambda - lambda.hat) ** 2))
      } else {
        low.ind.l <- ifelse(l == 1, 1, sum(n.vector[seq(l = 1)]) + 1)
        high.ind.l <- sum(n.vector[seq(l)])
        lambda.kl <- as.numeric(lambda.matrix[seq(low.ind.k, high.ind.k), l])
        lambda.lk <- as.numeric(lambda.matrix[seq(low.ind.l, high.ind.l), k])
        lambda.hats <- estimate.lambda.block(P.block)
        return(sum((lambda.kl - lambda.hats[[1]]) ** 2) + 
                 sum((lambda.lk - lambda.hats[[2]]) ** 2))
      }
    }) %>% 
      sum() %>% 
      return()
  }) %>% 
    sum() %>% 
    magrittr::divide_by(n * K) %>% 
    sqrt() %>% 
    return()
}