set.seed(475675)

source('~/dev/pabm-grdpg/functions.R')

n1 <- 200
n2 <- 200
n <- n1 + n2
z <- c(rep(1, n1), rep(2, n2))
p <- 1/50
q <- 1/50
r <- 1/5
omega <- rbeta(n, 2, 1)
Omega <- tcrossprod(omega)
P <- matrix(r, nrow = n, ncol = n)
P[seq(n1), seq(n1)] <- p
P[seq(n1 + 1, n), seq(n1 + 1, n)] <- q
P <- P * Omega
A <- draw.graph(P)
qgraph::qgraph(A, vsize = 4, groups = factor(z), legend = FALSE)

A.eigen <- eigen(A)
X.hat <- A.eigen$vectors[, c(1, n)] %*% diag(abs(A.eigen$values[c(1, n)]) ** .5)
plot(X.hat, 
     asp = 1,
     col = z * 2, 
     xlab = NA, ylab = NA)

norms <- apply(X.hat, 1, function(x) sqrt(sum(x ^ 2)))
cosine.sim <- X.hat %*% t(X.hat) / tcrossprod(norms)
eigen.cosine.sim <- eigen(cosine.sim, symmetric = TRUE)
Y <- sweep(eigen.cosine.sim$vectors[, 1:2], 2, 
           sqrt(eigen.cosine.sim$values[1:2]), `*`)
plot(Y, col = z * 2, asp = 1)
zhat <- kmeans(Y, 2)$cluster
table(z, zhat)

L <- diag(colSums(A)) - A
L.eigen <- eigen(L, symmetric = TRUE)
Y <- sweep(L.eigen$vectors[, c(n - 1, n - 2)], 2, sqrt(L.eigen$values[c(n - 1, n - 2)]), `/`)
plot(Y, asp = 1, col = z * 2)

em.out <- em.sbm(A, assortative = FALSE)
table(em.out$z, z)
