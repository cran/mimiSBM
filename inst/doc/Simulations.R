## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(mimiSBM)

## -----------------------------------------------------------------------------
set.seed(42)
N <- 100
V <- 10
K <- c(2,3,4)
pi_k <- rep(1/5,5) # Here the true number of clusters will be equal to 5
rho <- rep(1/length(K),length(K))
sorted = TRUE
p_switch = NULL

## -----------------------------------------------------------------------------
res <- rSMB_partition(N = N,V = V,K = K,pi_k = pi_k,rho = rho,sorted = sorted,p_switch = p_switch)

## -----------------------------------------------------------------------------
#plot_adjacency(res$simulation$A)
image(res$simulation$A[,,1])
image(res$simulation$A[,,3])
image(res$simulation$A[,,10])

## -----------------------------------------------------------------------------
(clust <- res$params$clusters)

## -----------------------------------------------------------------------------
length(clust)

## -----------------------------------------------------------------------------
clust[[1]]

## -----------------------------------------------------------------------------
sapply(1:3,function(s){image(res$params$alpha_klq[,,s])})

## -----------------------------------------------------------------------------
set.seed(42)
N <- 100
V <- 10
K <- c(2,3,4)
pi_k <- rep(1/5,5) # Here the true number of clusters will be equal to 5
rho <- rep(1/length(K),length(K))
sorted = TRUE
p_switch = 0.1

## -----------------------------------------------------------------------------
res <- rSMB_partition(N = N,V = V,K = K,pi_k = pi_k,rho = rho,sorted = sorted,p_switch = p_switch)

## -----------------------------------------------------------------------------
#plot_adjacency(res$simulation$A)
image(res$simulation$A[,,1])
image(res$simulation$A[,,3])
image(res$simulation$A[,,10])

## -----------------------------------------------------------------------------
res <- rSMB_partition(N = N,V = V,K = K,pi_k = pi_k,rho = rho,sorted = FALSE,p_switch = p_switch)
#plot_adjacency(res$simulation$A)
image(res$simulation$A[,,1])
image(res$simulation$A[,,3])
image(res$simulation$A[,,10])

