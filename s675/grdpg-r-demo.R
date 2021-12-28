data(karate, package = 'igraphdata')
A <- as.matrix(igraph::as_adjacency_matrix(karate, type = 'both'))
z <- igraph::vertex_attr(karate)$Faction
n <- length(z)
qgraph::qgraph(A)
eigen.A <- eigen(A, symmetric = TRUE)
V <- eigen.A$vectors
D <- eigen.A$values
X <- sweep(V[, c(1:3, n)], 2, sqrt(abs(D[c(1:3, n)])), `*`)
plot(X, asp = 1, col = z)
pairs(X, asp = 1, col = z)

V <- V[, c(1:3, n)]
B <- abs(V %*% t(V))

n <- 100
theta <- rbeta(n, 1/4, 1/4)
hist(theta)
X <- cbind(cos(pi * theta / 2 - pi / 4), sin(pi * theta / 2 - pi / 4))
plot(X, asp = 1)
P <- X %*% diag(c(1, -1)) %*% t(X)
summary(as.numeric(P))
A <- draw.graph(P)
qgraph::qgraph(A)
Xhat <- embedding(A, 1, 1)
plot(Xhat, asp = 1)
