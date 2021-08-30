embedding <- function(A, p = NULL, q = NULL,
                      scale = TRUE,
                      eps = 1e-6) {
  if (p + q == 0) {
    stop('one of p or q must be > 0')
  }
  n <- nrow(A)
  eigen.A <- eigen(A, symmetric = TRUE)
  if (is.null(p) | is.null(q)) {
    keep <- (abs(eigen.A$values) > eps)
  } else {
    if (p * q > 0) {
      keep <- c(seq(p), seq(n, n - q + 1))
    } else if (p == 0) {
      keep <- seq(n, n - q + 1)
    } else {
      keep <- seq(p)
    }
  }
  
  U <- eigen.A$vectors[, keep]
  if (scale) {
    S <- diag(sqrt(abs(eigen.A$values[keep])))
    return(U %*% S)
  } else {
    return(U * sqrt(n))
  }
}

draw.graph <- function(P) {
  n <- nrow(P)
  A <- matrix(0, nrow = n, ncol = n)
  A[upper.tri(A)] <- rbinom(n * (n - 1) / 2, 1, P[upper.tri(P)])
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
  B <- n * abs(V %*% t(V))
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

lambda.rmse <- function(P, A, clustering) {
  K <- max(clustering)
  n <- nrow(A)
  n.vector <- sapply(seq(K), function(k) sum(clustering == k))
  sapply(seq(K), function(k) {
    n.k <- sum(clustering == k)
    sapply(seq(k), function(l) {
      n.l <- sum(clustering == l)
      P.block <- P[clustering == k, clustering == l]
      A.block <- A[clustering == k, clustering == l]
      if (k == l) {
        # lambda <- estimate.lambda.block(P.block, TRUE)
        lambda.hat <- estimate.lambda.block(A.block, TRUE)
        P.hat <- lambda.hat %*% t(lambda.hat)
        # return(sum((lambda - lambda.hat) ** 2))
        return(norm(P.hat - P.block, 'F') ** 2)
      } else {
        # lambdas <- estimate.lambda.block(P.block)
        lambda.hats <- estimate.lambda.block(A.block)
        # return(sum((lambdas[[1]] - lambda.hats[[1]]) ** 2) + 
        #          sum((lambdas[[2]] - lambda.hats[[2]]) ** 2))
        P.hat <- lambda.hats[[1]] %*% t(lambda.hats[[2]])
        return(2 * norm(P.hat - P.block, 'F') ** 2)
      }
    }) %>% 
      sum() %>% 
      return()
  }) %>% 
    sum() %>% 
    # magrittr::divide_by(n * K) %>% 
    magrittr::divide_by(n ** 2) %>% 
    sqrt() %>%
    return()
}

ssc <- function(A, 
                K = 2,
                lambda = .01,
                parallel = FALSE,
                normalize = TRUE,
                scale = TRUE) {
  p <- K * (K + 1) / 2
  q <- K * (K - 1) / 2
  
  Y <- t(embedding(A, p, q, scale))
  
  if (normalize) {
    Y <- sweep(Y, 2, apply(Y, 2, function(y) sqrt(sum(y ** 2))), `/`)
    Y[is.nan(Y)] <- 0
  }
  
  N <- ncol(Y)
  B <- plyr::aaply(seq(N), 1, function(i) {
    y <- Y[, i]
    X <- Y[, -i]
    betahat <- glmnet::glmnet(X, y, lambda = lambda, intercept = FALSE) %>% 
      coef() %>% 
      as.numeric()
    if (i != N) {
      betahat <- c(betahat[seq(i)], 0, betahat[seq(i + 1, N)])
    } else {
      betahat <- c(betahat, 0)
    }
    betahat <- betahat[-1]
    return(betahat)
  }, .parallel = parallel) %>% 
    abs()
  B <- sweep(B, 2, apply(B, 2, max), `/`)
  B[is.nan(B)] <- 0
  W <- B + t(B)
  L <- normalized.laplacian(W)
  L.eigen <- eigen(L, symmetric = TRUE)
  X <- L.eigen$vectors[, seq(N, N - K)]
  # X <- sweep(X, 2, apply(X, 2, mean), `-`)
  # X <- sweep(X, 2, apply(X, 2, sd), `/`)
  clustering <- mclust::Mclust(X, K, verbose = FALSE)$classification
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
        # lambda <- estimate.lambda.block(P.kk, TRUE)
        # return(sum((lambda - lambda.hat) ** 2))
        P.hat <- lambda.hat %*% t(lambda.hat)
        return(norm(P.hat - P.kk, 'F') ** 2)
      } else {
        n.l <- n.vector[l]
        e.n.l <- rep(1, n.l)
        A.kl <- A[z == k, z == l]
        P.kl <- P[z == k, z == l]
        # lambdas <- estimate.lambda.block(P.kl)
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
        # return(sum((lambdas[[1]] - lambda.hat.kl) ** 2) + 
        #          sum((lambdas[[2]] - lambda.hat.lk) ** 2))
        P.hat <- lambda.hat.kl %*% t(lambda.hat.lk)
        return(2 * norm(P.kl - P.hat, 'F') ** 2)
      }
    }) %>% 
      sum() %>% 
      return()
  }) %>% 
    sum() %>% 
    # magrittr::divide_by(n * K) %>% 
    magrittr::divide_by(n ** 2) %>% 
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
    # theme_void() +
    # ggthemes::theme_base() + 
    theme(legend.position = 'none',
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    scale_x_continuous(expand = c(0, 0)) + 
    scale_y_continuous(expand = c(0, 0)) + 
    scale_fill_gradient(low = 'white', high = 'black')
  
  if (lines) {
    positions.x <- which(diff(z) != 0)
    positions.y <- which(diff(rev(z)) != 0)
    out.plot <- out.plot + 
      geom_hline(yintercept = positions.y + .5, 
                 colour = 'red',
                 size = .1) + 
      geom_vline(xintercept = positions.x + .5, 
                 colour = 'red',
                 size = .1)
  }
  
  return(out.plot)
}

mod.max.exchange <- function(A, init.clust) {
  K <- max(init.clust)
  n <- nrow(A)
  indices <- seq(n)
  clustering <- init.clust
  clusters <- seq(K)
  current.modularity <- Q.PA(A, clustering)
  prev.modularity <- current.modularity - 1
  print(current.modularity)
  
  iter <- 0
  
  while (prev.modularity != current.modularity) {
    prev.modularity <- current.modularity
    # indices <- sample(indices)
    for (i in indices) {
      current.label <- clustering[i]
      neighbors <- which(A[i, ] == 1)
      neighboring.labels <- unique(clustering[neighbors])
      if ((length(neighboring.labels) > 1) | 
          (neighboring.labels[1] != current.label)) {
        # print(i)
        for (k in neighboring.labels) {
          if (k != current.label) {
            proposed.clustering <- clustering
            proposed.clustering[i] <- k
            proposed.modularity <- Q.PA(A, proposed.clustering)
            if (proposed.modularity > current.modularity) {
              clustering <- proposed.clustering
              current.modularity <- proposed.modularity
            }
          }
        }
      }
    }
    iter <- iter + 1
    print(iter)
    print(current.modularity)
  }
  return(clustering)
}

ssc2 <- function(A, 
                 K = 2,
                 lambda = .01,
                 parallel = FALSE,
                 normalize = TRUE) {
  
  N <- ncol(A)
  B <- plyr::aaply(seq(N), 1, function(i) {
    y <- A[, i]
    X <- A[, -i]
    betahat <- glmnet::glmnet(X, y, lambda = lambda, intercept = FALSE) %>% 
      coef() %>% 
      as.numeric()
    if (i != N) {
      betahat <- c(betahat[seq(i)], 0, betahat[seq(i + 1, N)])
    } else {
      betahat <- c(betahat, 0)
    }
    betahat <- betahat[-1]
    return(betahat)
  }, .parallel = parallel) %>% 
    abs()
  if (normalize) {
    B <- sweep(B, 2, apply(B, 2, max), `/`)
    B[is.nan(B)] <- 0
  }
  W <- B + t(B)
  L <- normalized.laplacian(W)
  L.eigen <- eigen(L, symmetric = TRUE)
  X <- L.eigen$vectors[, seq(N, N - K)]
  clustering <- mclust::Mclust(X, K, verbose = FALSE)$classification
  return(clustering)
}
