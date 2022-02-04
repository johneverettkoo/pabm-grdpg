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

cluster.pabm <- function(A, K = 2, 
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

ssc <- function(A, 
                K = 2,
                lambda = .01,
                parallel = FALSE,
                normalize = TRUE,
                scale = TRUE, 
                d = K + 1,
                offset = 0) {
  # subspace clustering on the ASE of a PABM
  
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
  X <- L.eigen$vectors[, seq(N - offset, N - d + 1 - offset)]
  clustering <- mclust::Mclust(X, K, verbose = FALSE)$classification
  return(clustering)
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

mod.max <- function(A) {
  Ag <- igraph::graph_from_adjacency_matrix(A, mode = 'undirected')
  clustering <- igraph::cluster_louvain(Ag)
  return(igraph::membership(clustering))
}

cluster.acc <- function(yhat, yobs, reorder.mat) {
  K <- max(c(yhat, yobs))
  n <- length(yhat)
  if (missing(reorder.mat)) {
    reorder.mat <- gtools::permutations(K, K)
  }
  
  original <- seq(K)
  accuracies <- apply(reorder.mat, 1, function(reorder) {
    yhat.remap <- plyr::mapvalues(yhat, original, reorder)
    table(yhat.remap, yobs) %>% 
      as.matrix() %>% 
      diag() %>% 
      sum() %>% 
      magrittr::divide_by(n) %>% 
      return()
  })
  return(max(accuracies))
}

ssc2 <- function(A, 
                 K = 2,
                 lambda = .01,
                 parallel = FALSE,
                 normalize = TRUE) {
  # subspace clustering on the adjacency matrix
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

lambda.rmse <- function(P, A, clustering, rho = 1) {
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
    magrittr::divide_by(rho) %>% 
    return()
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

lambda.rmse.mle <- function(P, A, z, rho = 1) {
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
    magrittr::divide_by(rho) %>% 
    return()
}

cluster.sg <- function(A, K = 2) {
  # adapted from code provided by sengupta and chen:
  # a blockmodel for node popularity in networks with community structure
  
  f.PA<-function(A,b){	
    K<-max(b)       # no. of communities
    N<-nrow(A)      # no. of nodes
    M<-matrix(NA,nrow=N,ncol=K)  # popularity matrix
    O<-matrix(NA,nrow=K,ncol=K)  # community interaction matrix
    for (i in 1:N){		# calculate M
      for (r in 1:K){
        nodes = which(b == r)
        M[i,r] = sum(A[i,nodes])
      }}
    for (r in 1:K){		# calculate O
      for (s in r:K){
        nodes1 = which(b == r)
        nodes2 = which(b == s)
        O[r,s] = sum(A[nodes1,nodes2])
        O[s,r] = O[r,s]
      }}
    list(M=M, O=O)
  }
  
  Q.PA <- function(A, b){
    foo<-f.PA(A,b)
    O=foo$O; M = foo$M
    s1 = sum(M*log(M),na.rm=TRUE) # na.rm = TRUE ignores M=0 cases as log(0) = NA
    s2 = sum(O*log(O),na.rm=TRUE) # na.rm = TRUE ignores O=0 cases as log(0) = NA
    return(2*s1-s2)
  }
  
  f.DC<-function(A,b){	
    K<-max(b)
    O<-matrix(NA,nrow=K,ncol=K)  # interaction matrix
    for (i in 1:K){		# calculate O
      for (j in i:K){
        nodes1 = which(b == i)
        nodes2 = which(b == j)
        O[i,j] = sum(A[nodes1,nodes2])
        O[j,i] = O[i,j]
      }}
    list(O=O, d=rowSums(O))
  }
  
  Q.DC <- function(A, b){
    K <- max(b)
    q <- matrix(0, nrow = K, ncol = K)
    foo<-f.DC(A,b)
    O<-foo$O; d<-foo$d 
    for (i in 1:K){
      for (j in 1:K){
        if (O[i,j]>0){
          q[i,j] = O[i,j]*log(O[i,j]/(d[i]*d[j]))}
      }}
    return(sum(q))
  }	# formula for Q
  
  b.can = EPalgo(A,eps=0, K = K) # EP algorithm (no perturbation)
  Q.PA.can = rep(NA, ncol(b.can))	# array to store Q values
  Q.DC.can = rep(NA, ncol(b.can))	# array to store Q values
  for (i in 1:ncol(b.can)){
    #check if any cluster is empty
    foo = rep(NA, K)
    for (clus in 1:K) {foo[clus]=sum(b.can[,i]==clus)}
    if (min(foo)==0) {stop('Empty groups are not allowed')} 
    Q.PA.can[i] = Q.PA(A, b=b.can[,i])   # fit PABM
    Q.DC.can[i] = Q.DC(A, b=b.can[,i])   # fit DCBM
  } # end of i for loop
  foo1 = order(-Q.PA.can)[1] 
  b.PA = b.can[,foo1]   # community assignment that maximises Q.PA
  return(b.PA)
}

EPalgo<-function(A,eps=0, K = 2){
  # provided by sengupta and chen:
  # a blockmodel for node popularity in networks with community structure
  
  ##### perturbed adj matrix ##### (Le pg 15, Amini 2013)
  tau = eps*(mean(colSums(A))/nrow(A))
  A = A + tau
  
  foo<-eigen(A, symmetric = TRUE)
  val = abs(foo$values)			# pick the top 2
  id = order(-val)				# eigenvalues of A and put
  id_vec = id[1:K]				# their eigenvectors into a 2*N matrix
  # columns of foo$vectors = eigenvectors of A
  # we want a 2xn matrix whose rows are the leading eigenvectors
  X = t(foo$vectors[,id_vec])		
  y = X[,1:K]
  
  comms = 1:K
  u <- list(comms)
  v = expand.grid(rep(u, K))
  v = as.matrix(v)
  # initialize with the parallelogram
  epts = y%*%t(v)	# extreme pts are the columns of this matrix
  b.can = t(v)	# candidate configurations.
  row.names(b.can)=NULL
  
  ptm<-proc.time()
  for (i in (K+1):ncol(X)){
    # b.can1 = rbind(b.can,rep(1,ncol(b.can)))
    # b.can2 = rbind(b.can,rep(2,ncol(b.can)))
    # b.can = cbind(b.can1,b.can2)
    b.can.list <- lapply(seq(K), function(k) rbind(b.can, rep(k, ncol(b.can))))
    b.can <- do.call(cbind, b.can.list)
    foo = X[,1:i]%*%b.can
    hull = chull(t(foo))
    epts = foo[,hull]
    b.can = b.can[,hull]
  }	# next i = next row of X
  proc.time()-ptm
  
  ##### remove invalid candidates
  k = max(b.can)
  foo = b.can
  foo1 = NA
  for (i in 1:ncol(b.can)){
    foo2 = rep(NA,k)
    for (clus in 1:k){foo2[clus]=sum(b.can[,i]==clus)}
    if (min(foo2)==0){foo1 = c(foo1,i)}
  }
  if (length(foo1)>1) {foo1 = foo1[-1]
  b.can = b.can[,-foo1]}
  
  ###### remove eqv candidates
  foo1 = NA
  for (i in 2:ncol(b.can)){
    for (j in 1:i){
      foo4 = abs(b.can[,i] - b.can[,j])
      if (mean(foo4) == 1){ # this means b.can[,i] and b.can[,j] are exactly 1 apart
        foo1 = c(foo1,i)
        break}
    }}
  if (length(foo1)>1){
    foo1 = foo1[-1]
    b.can = b.can[,-foo1]
  }
  return(b.can)
}	# end of function
