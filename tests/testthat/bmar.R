library(huge)
devtools::load_all("~/bmar")

sim <- huge.generator(n = 20, d = 10, v = 0.3, u = 0.1, graph = "random")

x <- sim$data
l <- 250
niter <- 50
burnin <- 0

bmar(x = x, l = l, niter = niter, burnin = burnin)
