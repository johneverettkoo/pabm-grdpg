import::from(magrittr, `%>%`)
import::from(readr, read_csv)
library(mclust)
library(haven)
library(ggplot2)

type <- 'visitgo'
villageno <- 0

# villageno <- villageno + 1

classification <- 'hohreligion'
keep <- c(1, 2)
villageno <- 12

setwd('~/dev/multiresolution_networks/data/indian_village_raw/datav4.0/Data')

villageno <- villageno + 1
print(villageno)

A <- read_csv(paste0("1. Network Data/Adjacency Matrices/adj_", type, "_HH_vilno_", villageno, ".csv"), 
              col_names = FALSE) %>% 
  as.matrix() %>% 
  magrittr::set_colnames(NULL)
deg <- colSums(A)
A <- A[deg > 0, deg > 0]
household_characteristics <- haven::read_dta("2. Demographics and Outcomes/household_characteristics.dta") %>% 
  dplyr::filter(village == villageno)
household_characteristics <- household_characteristics[deg > 0, ]
z <- as.numeric(household_characteristics[[classification]])
A <- A[z %in% keep, z %in% keep]
z <- z[z %in% keep]
deg <- colSums(A)
A <- A[deg > 0, deg > 0]
z <- z[deg > 0]
table(z)
qgraph::qgraph(A, groups = factor(z))

saveRDS(list(A = A, z = z), 
        paste0('~/dev/pabm-grdpg/data/village-hh-', villageno, '.rds'))

K <- length(unique(z))

p <- 3
q <- 1
normalize <- FALSE
d.eigenmap <- K + 1
laplacian <- 'normalized'
use.all <- TRUE
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
  eigenmap <- eigen(L, symmetric = TRUE)$vectors[, seq(n, n - d.eigenmap + 1)]
} else {
  eigenmap <- eigen(L, symmetric = TRUE)$vectors[, seq(n - 1, n - d.eigenmap + 1)]
}
pairs(eigenmap, col = factor(z))
zhat <- mclust::Mclust(eigenmap, K)$classification
table(z, zhat)
plot.A(A, zhat)
plot.A(A, z)
qgraph::qgraph(A, groups = factor(zhat))

zhat.ssc <- ssc(A, lambda = .1)
table(zhat.ssc, z)
