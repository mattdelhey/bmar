library(huge)
library(xtable)
library(parallel)

cl <- makeCluster(4)
clusterEvalQ(cl, {library(huge); library(xtable); devtools::load_all("~/bmar")})

## Parameters
f.out.graphs <- "simulation-graphs.pdf"
reps <- 3

## Simulations
sim.plot(sim.objects, f.out.graphs)

## sim.result[method, error, model, replicate]
## sim.results <- simplify2array(replicate(reps, {
##     sim.random  <- huge.generator(n = 200, d = 300, v = 0.3, u = 0.1, graph = "random", prob = 0.006)
##     sim.cluster <- huge.generator(n = 200, d = 300, v = 0.3, u = 0.1, graph = "cluster")
##     sim.band    <- huge.generator(n = 200, d = 300, v = 0.3, u = 0.1, graph = "band")
##     sim.sclfree <- huge.generator(n = 200, d = 300, graph = "scale-free")    
##     sim.objects <- list(sim.random, sim.cluster, sim.band, sim.sclfree)
##     sim.summarize(sim.objects, nlambda = 100, K = 5, stars.thresh = 0.05, stars.subsample.ratio = NULL, rep.num = 10)
## }))

sim.results <- simplify2array(parLapply(cl, 1:reps, function(i) {
    sim.random  <- huge.generator(n = 200, d = 300, v = 0.3, u = 0.1, graph = "random", prob = 0.006)
    sim.cluster <- huge.generator(n = 200, d = 300, v = 0.3, u = 0.1, graph = "cluster")
    sim.band    <- huge.generator(n = 200, d = 300, v = 0.3, u = 0.1, graph = "band")
    sim.sclfree <- huge.generator(n = 200, d = 300, graph = "scale-free")    
    sim.objects <- list(sim.random, sim.cluster, sim.band, sim.sclfree)
    sim.summarize(cl, sim.objects, nlambda = 30, K = 5, stars.thresh = 0.05, stars.subsample.ratio = NULL, rep.num = 10)
}))

## Mean, sd
means <- apply(sim.results, c(1, 2, 3), mean, na.rm = TRUE)
apply(means[, c(-3, -4), ], 3, xtable)

sds <- apply(sim.results, c(1, 2, 3), sd, na.rm = TRUE)
apply(sds[, c(-3, -4), ], 3, xtable)
