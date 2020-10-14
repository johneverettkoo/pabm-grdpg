embedding <- function(A, p = NULL, q = NULL,
                      eps = 1e-6) {
  n <- nrow(A)
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

graph.laplacian <- function(W) {
  D <- diag(colSums(W))
  L <- D - W
  return(L)
}

cluster.pabm <- function(A, K, 
                         normalize = TRUE, 
                         use.all = TRUE,
                         p = NULL, q = NULL, d.eigenmap = K + 1, 
                         laplacian = 'normalized') {
  n <- nrow(A)
  if (is.null(p)) {
    p <- K * (K + 1) / 2
  }
  if (is.null(q)) {
    q <- K * (K - 1) / 2
  }
  if (p * q > 0) {
    indices <- c(seq(p), seq(n, n - q + 1))
  } else if (p == 0) {
    indices <- seq(n, n - q + q)
  } else if (q == 0) {
    indices <- seq(p)
  }
  V <- eigen(A, symmetric = TRUE)$vectors[, indices]
  if (normalize) {
    v.norm <- sqrt(rowSums(V ** 2))
    V <- sweep(V, 1, v.norm, `/`)
  }
  B <- abs(V %*% t(V))
  if (laplacian == 'graph') {
    L <- graph.laplacian(B)
  } else if (laplacian == 'normalized') {
    L <- normalized.laplacian(B)
  }
  if (use.all) {
    eigenmap <- 
      eigen(L, symmetric = TRUE)$vectors[, seq(n, n - d.eigenmap + 1)]
  } else {
    eigenmap <- 
      eigen(L, symmetric = TRUE)$vectors[, seq(n - 1, n - d.eigenmap + 1)]
  }
  clustering <- mclust::Mclust(eigenmap, K)$classification
  
  return(clustering)
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

lambda.rmse <- function(P, Phat, clustering) {
  K <- max(clustering)
  n <- nrow(Phat)
  n.vector <- sapply(seq(K), function(k) sum(clustering == k))
  sapply(seq(K), function(k) {
    sapply(seq(k), function(l) {
      P.block <- P[clustering == k, clustering == l]
      Phat.block <- Phat[clustering == k, clustering == l]
      if (k == l) {
        lambda <- estimate.lambda.block(P.block, TRUE)
        lambda.hat <- estimate.lambda.block(Phat.block, TRUE)
        return(sum((lambda - lambda.hat) ** 2))
      } else {
        lambdas <- estimate.lambda.block(P.block)
        lambda.hats <- estimate.lambda.block(Phat.block)
        return(sum((lambdas[[1]] - lambda.hats[[1]]) ** 2) + 
                 sum((lambdas[[2]] - lambda.hats[[2]]) ** 2))
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

ssc <- function(A, 
                K = 2,
                lambda = .01,
                parallel = FALSE) {
  p <- K * (K + 1) / 2
  q <- K * (K - 1) / 2
  
  Y <- t(embedding(A, p, q))
  
  n <- nrow(Y)
  N <- ncol(Y)
  B <- plyr::aaply(seq(N), 1, function(i) {
    y <- Y[, i]
    X <- Y[, -i]
    betahat <- glmnet::glmnet(X, y, lambda = lambda) %>% 
      coef() %>% 
      as.numeric()
    return(betahat)
  }, .parallel = parallel) %>% 
    abs()
  W <- B + t(B)
  L <- normalized.laplacian(W)
  L.eigen <- eigen(L, symmetric = TRUE)
  X <- diag(colSums(W) ** -.5) %*% L.eigen$vectors[, seq(N - 1, N - K)]
  clustering <- mclust::Mclust(X, K)$classification
  return(clustering)
}

lambda.rmse.mle <- function(P, A, z) {
  K <- max(z)
  n <- nrow(A)
  n.vector <- sapply(seq(K), function(k) sum(z == k))
  sapply(seq(K), function(k) {
    n.k <- n.vector[k]
    e.n.k <- rep(1, n.k)
    sapply(seq(k), function(l) {
      if (k == l) {
        A.kk <- A[z == k, z == l]
        P.kk <- P[z == k, z == l]
        lambda.hat <- 
          as.numeric((A.kk %*% e.n.k) / 
                       as.numeric(sqrt(t(e.n.k) %*% A.kk %*% e.n.k)))
        lambda <- estimate.lambda.block(P.kk, TRUE)
        return(sum((lambda - lambda.hat) ** 2))
      } else {
        n.l <- n.vector[l]
        e.n.l <- rep(1, n.l)
        A.kl <- A[z == k, z == l]
        P.kl <- P[z == k, z == l]
        lambdas <- estimate.lambda.block(P.kl)
        lambda.hat.kl <- 
          as.numeric((A.kl %*% e.n.l) / 
                       as.numeric(sqrt(t(e.n.k) %*% A.kl %*% e.n.l)))
        lambda.hat.lk <- 
          as.numeric((t(A.kl) %*% e.n.k) / 
                       as.numeric(sqrt(t(e.n.k) %*% A.kl %*% e.n.l)))
        lambda.hat.kl <- dplyr::if_else(is.nan(lambda.hat.kl),
                                        0, lambda.hat.kl)
        lambda.hat.lk <- dplyr::if_else(is.nan(lambda.hat.lk), 
                                        0, lambda.hat.lk)
        return(sum((lambdas[[1]] - lambda.hat.kl) ** 2) + 
                 sum((lambdas[[2]] - lambda.hat.lk) ** 2))
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

generate.P.beta <- function(n, K = 2, a1 = 2, b1 = 1, a2 = 1, b2 = 2,
                            unbalanced = FALSE) {
  if (unbalanced) {
    clustering <- rmultinom(n, 1, seq(K) ** -1) %>% 
      apply(2, function(x) which(x == 1))
  } else {
    clustering <- sample(seq(K), n, replace = TRUE)
  }
  clustering <- sort(clustering)
  n.vector <- sapply(seq(K), function(k) sum(clustering == k))
  
  P <- matrix(NA, n, n)
  for (k in seq(K)) {
    n.k <- n.vector[k]
    low.ind.k <- ifelse(k == 1, 1, sum(n.vector[seq(k - 1)]) + 1)
    high.ind.k <- sum(n.vector[seq(k)])
    for (l in seq(k)) {
      n.l <- n.vector[l]
      if (k == l) {
        lambda <- rbeta(n.k, a1, b1)
        P.kk <- lambda %*% t(lambda)
        P[low.ind.k:high.ind.k, low.ind.k:high.ind.k] <- P.kk
      } else {
        low.ind.l <- ifelse(l == 1, 1, sum(n.vector[seq(l - 1)]) + 1)
        high.ind.l <- sum(n.vector[seq(l)])
        lambda.kl <- rbeta(n.k, a2, b2)
        lambda.lk <- rbeta(n.l, a2, b2)
        P.kl <- lambda.kl %*% t(lambda.lk)
        P.lk <- lambda.lk %*% t(lambda.kl)
        P[low.ind.k:high.ind.k, low.ind.l:high.ind.l] <- P.kl
        P[low.ind.l:high.ind.l, low.ind.k:high.ind.k] <- P.lk
      }
    }
  }
  return(list(P = P, clustering = clustering))
}

plot.A <- function(A, z, lines = TRUE, max.size = 500) {
  # https://www.r-graph-gallery.com/79-levelplot-with-ggplot2.html
  A <- A[order(z), order(z)]
  z <- sort(z)
  n <- nrow(A)
  if (n > max.size) {
    ind <- sort(sample(seq(n), max.size))
    A <- A[ind, ind]
    z <- z[ind]
  }
  out.plot <- A %>% 
    tibble::as_tibble() %>% 
    tibble::rowid_to_column(var = 'X') %>% 
    tidyr::gather(key = 'Y', value = 'Z', -1) %>% 
    dplyr::mutate(Y = rev(as.numeric(gsub('V', '', Y)))) %>% 
    ggplot() + 
    geom_tile(aes(x = X, y = Y, fill = Z)) + 
    coord_fixed() + 
    labs(x = NULL, y = NULL) + 
    theme_void() +
    theme(legend.position = 'none',
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    scale_x_continuous(expand = c(0, 0)) + 
    scale_y_continuous(expand = c(0, 0)) + 
    scale_fill_gradient(low = 'red', high = 'yellow')
  
  if (lines) {
    positions.x <- which(diff(z) != 0)
    positions.y <- which(diff(rev(z)) != 0)
    out.plot <- out.plot + 
      geom_hline(yintercept = positions.y + .5) + 
      geom_vline(xintercept = positions.x + .5)
  }
  
  return(out.plot)
}
