# karanataka villages network data
# https://github.com/tedwestling/multiresolution_networks.git

library(mclust)
import::from(magrittr, `%>%`)
setwd('~/dev/pabm-grdpg')
source('simulations-and-data-analysis/functions.R')

remove.deg.0 <- function(A, z, keep = unique(z)) {
  T <- colSums(A)
  A <- A[T > 0, T > 0]
  z <- z[T > 0]
  A <- A[z %in% keep, z %in% keep]
  z <- z[z %in% keep]
  return(list(A = A, z = z))
}

compute.errors <- function(A, z, K = 2, lambda = 1e-1) {
  zhat.osc <- cluster.pabm(A, K)
  zhat.ssc <- ssc(A, K, lambda)
  zhat.mm <- cluster.sg(A)
  list(error.osc = 1 - cluster.acc(z, zhat.osc),
       error.ssc = 1 - cluster.acc(z, zhat.ssc),
       error.mm = 1 - cluster.acc(z, zhat.mm))
}

adj.path <- 'datav4.0/Data/1. Network Data/Adjacency Matrices/'
labels.df <- haven::read_dta(file.path('datav4.0/Data', 
                                       '2. Demographics and Outcomes',
                                       'household_characteristics.dta'))

# village 12
A12 <- readr::read_csv(file.path(adj.path, 'adj_visitgo_HH_vilno_12.csv'),
                       col_names = FALSE) %>% 
  as.matrix()
z12 <- labels.df %>% 
  dplyr::filter(village == 12) %>%
  .$hohreligion %>% 
  as.numeric()
Az12 <- remove.deg.0(A12, z12, c(1, 2))
A12 <- Az12$A
z12 <- Az12$z
compute.errors(A12, z12)

# village 31
A31 <- readr::read_csv(file.path(adj.path, 'adj_visitgo_HH_vilno_31.csv'),
                       col_names = FALSE) %>% 
  as.matrix()
z31 <- labels.df %>% 
  dplyr::filter(village == 31) %>%
  .$hohreligion %>% 
  as.numeric()
Az31 <- remove.deg.0(A31, z31, c(1, 2))
A31 <- Az31$A
z31 <- Az31$z
compute.errors(A31, z31)

# village 46
A46 <- readr::read_csv(file.path(adj.path, 'adj_visitgo_HH_vilno_46.csv'),
                       col_names = FALSE) %>% 
  as.matrix()
z46 <- labels.df %>% 
  dplyr::filter(village == 46) %>%
  .$hohreligion %>% 
  as.numeric()
Az46 <- remove.deg.0(A46, z46, c(1, 2))
A46 <- Az46$A
z46 <- Az46$z
compute.errors(A46, z46)
