import::from(magrittr, `%>%`)
library(plot.matrix)
import::from(foreach, `%do%`)
# library(CVXR)

source('http://pages.iu.edu/~mtrosset/Courses/675/stress.r')

doMC::registerDoMC(parallel::detectCores())

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

n1 <- 2 ** 5
n2 <- n1
n <- n1 + n2

z <- c(rep(1, n1), 
       rep(2, n2))

Ipq <- diag(c(1, 1, 1, -1))

lambda11 <- rbeta(n1, 2, 1)
lambda22 <- rbeta(n2, 2, 1)
lambda12 <- rbeta(n1, 1, 2)
lambda21 <- rbeta(n2, 1, 2)

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

Z <- embedding(P)
Zhat <- embedding(A, 3, 1)
pairs(Z, col = z, asp = 1)

P - Z %*% Ipq %*% t(Z)

XU <- cbind(c(lambda11, rep(0, n2)),
            c(rep(0, n1), lambda22),
            c(lambda12 / sqrt(2), lambda21 / sqrt(2)),
            c(-lambda12 / sqrt(2), lambda21 / sqrt(2)))
pairs(XU, col = z, asp = 1)
summary(as.numeric(P - XU %*% Ipq %*% t(XU)))
summary(as.numeric(P - Z %*% Ipq %*% t(Z)))

pairs(X, col = z, asp = 1)

perm <- cbind(c(1, 0, 0, 0),
              c(0, 0, 1, 0),
              c(0, 1, 0, 0),
              c(0, 0, 0, 1))

X %*% perm %*% t(X)
Z %*% perm %*% t(Z)

U <- eigen(perm)$vectors[, c(3, 1, 2, 4)]
Xhat <- Z %*% t(U)
pairs(Xhat, col = z, asp = 1)

summary(as.numeric(Xhat %*% U %*% Ipq %*% t(Xhat %*% U) - P))

normalized.laplacian <- function(W) {
  D <- diag(colSums(W) ** -.5)
  I <- diag(nrow(W))
  return(I - D %*% W %*% D)
}

ssc <- function(Y, 
                K = NULL,
                parallel = TRUE,
                lambda = .01) {
  
  n <- nrow(Y)
  N <- ncol(Y)
  B <- plyr::aaply(seq(N), 1, function(i) {
    y <- Y[, i]
    X <- Y[, -i]
    if (is.null(lambda)) {
      cv.lasso <- glmnet::cv.glmnet(X, y, grouped = FALSE)
      lambda <- cv.lasso$lambda.min
    }
    betahat <- glmnet::glmnet(X, y, lambda = lambda) %>% 
      coef() %>% 
      as.numeric()
    return(betahat)
  }, .parallel = parallel) %>% 
    abs()
  W <- B + t(B)
  L <- normalized.laplacian(W)
  L.eigen <- eigen(L, symmetric = TRUE)
  if (is.null(K)) {
    warning('unknown number of clusters--estimating ...')
    K <- N - which.max(-diff(deltas))
  }
  
  X <- diag(colSums(W) ** -.5) %*% L.eigen$vectors[, seq(N - 1, N - K)]
  # clustering <- kmeans(X, K)$cluster
  clustering <- mclust::Mclust(X, K)$classification
  return(clustering)
  
  # out <- foreach::foreach(j = seq(K), .combine = cbind) %do% {
  #   ind <- which(clustering == j)
  #   x <- t(Y)[ind, ]
  #   x.mean <- colMeans(x)
  #   x.tilde <- sweep(x, 2, x.mean)
  #   x.eigen <- eigen(x %*% t(x), symmetric = TRUE)
  #   x.proj <- sweep(x.eigen$vectors[, seq(K)], 
  #                   2,
  #                   x.eigen$values[seq(K)] ** .5,
  #                   `/`)
  #   if (j == 1) {
  #     rbind(x.proj, matrix(0, nrow = sum(clustering != j), ncol = K))
  #   } else if (j == K) {
  #     rbind(matrix(0, nrow = sum(clustering != j), ncol = K),
  #           x.proj)
  #   } else {
  #     rbind(matrix(0, nrow = sum(clustering < j), ncol = K),
  #           x.proj,
  #           matrix(0, nrow = sum(clustering > j), ncol = K))
  #   }
  # }
  # return(out)
}

pabm.clust.step <- function(A, K, eps = 1e-3) {
  p <- K * (K + 1) / 2
  q <- K * (K - 1) / 2
  n <- nrow(A)
  
  eigen.A <- eigen(A, symmetric = TRUE)
  v <- eigen.A$vectors[, c(seq(p), seq(n, n - q + 1))]
  W <- ((abs(v %*% t(v))) > eps) * 1
  diag(W) <- 0
  L <- graph.laplacian(W)
  eigen.L <- eigen(L, symmetric = TRUE)
}

car::scatter3d(Z[, 1], Z[, 2], Z[, 3],
               surface = FALSE, 
               point.col = z)

B <- lsa::cosine(t(Z))
plot(density(B))
hist(B)

B <- lsa::cosine(t(XU))

Q.inv <- MASS::ginv(XU) %*% Z
Q <- solve(Q.inv)

out <- vegan::procrustes(XU, Z)

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

min.vec <- function(x, a) sapply(x, function(y) min(y, a))

obj.fun <- function(angles, Z.hat) {
  R <- construct.Q(angles)
  # rot.angles <- angles[1:3]
  # hyp.angles <- angles[4:6]
  # hyp.angles <- c(angles[1], angles[4], -angles[4])
  # rot.angles <- c(angles[1], angles[2], -angles[2])
  # hyp.angles <- c(angles[3], angles[4], -angles[4])
  # rot.matrices <- lapply(seq(3), function(i) rot.mat.2(rot.angles[i], i))
  # hyp.matrices <- lapply(seq(3), function(i) hyp.rot.mat.2(hyp.angles[i], i))
  # matrices <- c(rot.matrices, hyp.matrices)
  # R <- Reduce(`%*%`, matrices)
  
  XU.c <- Z.hat %*% R
  
  # return(
  #   sum(apply(XU.c[, 1:2] ** 2, 1, min)) +
  #     sum(apply(cbind((XU.c[, 3] - XU.c[, 4]) ** 2,
  #                     (XU.c[, 3] + XU.c[, 4]) ** 2),
  #               1, min)) 
  # )
  clusters <- apply(XU.c[, 1:2] ** 2, 1, which.min)
  in.prods <- XU.c[clusters == 1, ] %*% t(XU.c[clusters == 2, , drop = FALSE])
  # return(
  #   sum(apply(XU.c[, 1:2] ** 2, 1, min)) +
  #     sum(apply(cbind((XU.c[, 3] - XU.c[, 4]) ** 2,
  #                     (XU.c[, 3] + XU.c[, 4]) ** 2),
  #               1, min)) + 
  #     Matrix::norm(in.prods, 'F')
  #   )
  # return(
  #   sum(apply(XU.c[, 1:2] ** 2, 1, min)) + 
  #     sum((XU.c[apply(XU.c[, 1:2] ** 2, 1, which.min) == 1, 3] + 
  #            XU.c[apply(XU.c[, 1:2] ** 2, 1, which.min) == 1, 4]) ** 2) * 
  #     sum(abs(XU.c[apply(XU.c[, 1:2] ** 2, 1, which.min) == 1, 3])) + 
  #     sum((XU.c[apply(XU.c[, 1:2] ** 2, 1, which.min) == 2, 3] - 
  #            XU.c[apply(XU.c[, 1:2] ** 2, 1, which.min) == 2, 4]) ** 2) * 
  #     sum(abs(XU.c[apply(XU.c[, 1:2] ** 2, 1, which.min) == 1, 3]))
  # )
  return(
    sum(apply(XU.c[, 1:2] ** 2, 1, min)) +
      sum((XU.c[clusters == 1, 3] +
             XU.c[clusters == 1, 4]) ** 2) +
      sum((XU.c[clusters == 2, 3] -
             XU.c[clusters == 2, 4]) ** 2) +
      Matrix::norm(in.prods,
                   'F')
  )
}

obj.fun(rep(0, 4), XU)

angles <- rep(0, 4)
out <- optim(angles, function(x) obj.fun(x, Z), method = 'Nelder-Mead',
             control = list(maxit = 1e5, abstol = 1e-10))
out <- optim(angles, function(x) obj.fun(x, Z), method = 'SANN',
             control = list(maxit = 1e5))
out <- optim(rep(0, 6), function(x) obj.fun(x, Z), method = 'L-BFGS-B',
             lower = 0, upper = 2 * pi, 
             control = list(maxit = 1e5))

angles <- rep(0, 4)
out <- optim(angles, function(x) obj.fun(x, Zhat), method = 'SANN',
             control = list(maxit = 1e5))
out <- optim(out$par, function(x) obj.fun(x, Zhat), method = 'Nelder-Mead',
             control = list(maxit = 1e5, abstol = 1e-10))

optim(runif(4, -pi, pi), function(x) {
  rot.angles <- x[1:3]
  # rot.angles <- x[4:6]
  hyp.angles <- c(x[1], x[4], -x[4])
  # rot.angles <- c(x[1], -x[2], x[2])
  # hyp.angles <- c(x[3], x[4], -x[4])
  rot.matrices <- lapply(seq(3), function(i) rot.mat.2(rot.angles[i], i))
  hyp.matrices <- lapply(seq(3), function(i) hyp.rot.mat.2(hyp.angles[i], i))
  matrices <- c(rot.matrices, hyp.matrices)
  R <- Reduce(`%*%`, rev(matrices))
  Matrix::norm(R - Q, 'F')
}, 
# control = list(abstol = 1e-10),
method = 'L-BFGS-B')

obj.fun <- function(angles, Z.hat, z, eps = 1e-6) {
  rot.angles <- angles[1:3]
  # hyp.angles <- angles[4:6]
  hyp.angles <- c(angles[1], angles[4], -angles[4])
  rot.matrices <- lapply(seq(3), function(i) rot.mat.2(rot.angles[i], i))
  hyp.matrices <- lapply(seq(3), function(i) hyp.rot.mat.2(hyp.angles[i], i))
  matrices <- c(rot.matrices, hyp.matrices)
  R <- Reduce(`%*%`, matrices)
  
  XU.c <- Z.hat %*% R
  
  return(
    sum(XU.c[z == 1, 1] ** 2) + 
      sum(XU.c[z == 2, 2] ** 2) + 
      sum((XU.c[z == 1, 3] - XU.c[z == 1, 4]) ** 2) * 
      sum(XU.c[z == 2, 3] ** 2) + 
      sum((XU.c[z == 2, 3] + XU.c[z == 2, 4]) ** 2) * 
      sum(XU.c[z == 1, 3] ** 2)
  )
}

out <- optim(angles, function(x) obj.fun(x, Z, z), method = 'Nelder-Mead',
             control = list(maxit = 1e5))
out <- optim(rep(0, 6), function(x) obj.fun(x, Z, z), method = 'L-BFGS-B',
             lower = 0, upper = 2 * pi, 
             control = list(maxit = 1e5))
out <- optim(angles, function(x) obj.fun(x, Z, z), method = 'SANN',
             control = list(maxit = 1e5))

lorentz <- function(v1, v2, v3) {
  v.sq <- v1 ** 2 + v2 ** 2 + v3 ** 2
  
  if (v.sq == 0) {
    return(diag(c(1, 1, 1, -1)))
  }
  
  gam <- 1 / sqrt(1 - v.sq)
  g1 <- gam - 1
  return(matrix(c(
    1 + g1 * v1 ** 2 / v.sq, g1 * v1 * v2 / v.sq, g1 * v1 * v3 / v.sq, -gam * v1,
    g1 * v2 * v1 / v.sq, 1 + g1 * v2 ** 2 / v.sq, g1 * v2 * v3 / v.sq, -gam * v2,
    g1 * v3 * v1 / v.sq, g1 * v3 * v2 / v.sq, 1 + g1 * v3 ** 2 / v.sq, -gam * v3,
    -gam * v1, -gam * v2, -gam * v3, gam
  ), byrow = TRUE, nrow = 4, ncol = 4))
}

rot.3d <- function(angles) {
  # https://en.wikipedia.org/wiki/Rotation_matrix#In_three_dimensions
  a1 <- angles[1]
  a2 <- angles[2]
  a3 <- angles[3]
  matrix(c(
    cos(a3) * cos(a2), cos(a3) * sin(a2) * sin(a1) - sin(a3) * cos(a1), cos(a3) * sin(a2) * cos(a1) + sin(a3) * sin(a1),
    sin(a3) * cos(a2), sin(a3) * sin(a2) * sin(a1) + cos(a3) * cos(a1), sin(a3) * sin(a2) * cos(a1) - cos(a3) * sin(a1),
    -sin(a2), cos(a2) * sin(a1), cos(a2) * cos(a1)
  ), byrow = TRUE, nrow = 3, ncol = 3)
}

hom.lorentz <- function(angles, b, s = -1) {
  M <- rot.3d(angles)
  gam <- sqrt(as.numeric(t(b) %*% b) + 1) * s
  a <- M %*% b / gam
  
  R <- matrix(nrow = 4, ncol = 4)
  R[1:3, 1:3] <- M
  R[4, 4] <- gam
  R[1:3, 4] <- -a
  R[4, 1:3] <- -t(b)
  return(R)
}

hom.lorentz <- function(angles, gam, a) {
  M <- rot.3d(angles)
  
  b <- gam * solve(M) %*% a
  
  R <- matrix(nrow = 4, ncol = 4)
  R[1:3, 1:3] <- M
  R[4, 4] <- gam
  R[1:3, 4] <- -a
  R[4, 1:3] <- -t(b)
  return(R)
}

angle <- runif(1, -pi, pi)
S <- matrix(c(cosh(angle), sinh(angle), 0, 0,
              sinh(angle), cosh(angle), 0, 0,
              0, 0, 1, 0,
              0, 0, 0, 1), 
            byrow = TRUE, 
            nrow = 4, ncol = 4)
