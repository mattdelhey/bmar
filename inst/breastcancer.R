library(huge)
devtools::install_github("mattdelhey/bmar")
library(bmar)

data(breastcancer, package = "gRbase")
bc <- as.matrix(breastcancer[, -1001]) # remove label

g <- huge(bc, nlambda = 100, method = "glasso", lambda.min.ratio = 0.01)
g2 <- huge(bc, nlambda = 100, method = "glasso", lambda.min.ratio = 0.01)
g3 <- huge(bc, nlambda = 10, method = "glasso", lambda.min.ratio = 0.01)
lambda.seq <- g$lambda

g.bic <- select.bic(g = g, n = nrow(bc), d = ncol(bc))
g.ss <- huge.select(g, criterion = "stars", stars.thresh = 0.05, stars.subsample.ratio = 0.5, rep.num = 5)
g.cv <- g.cvf(K = 5, x = bc, lambda.seq = lambda.seq)

g.bma <- bmar(g=g2, niter = 1000, l = 5, burnin = 200)
bma.edges <- and.rule(g.bma$theta.hat)

g.cv.min <- ose.rule(g.cv$m, g.cv$sd)

data.plot <- function(g = NULL, edges = NULL, index, main) {
    if (!is.null(g)) {
        g <- graph.adjacency(as.matrix(g$path[[index]]), mode = "undirected", diag = FALSE)
    }
    if (!is.null(edges)) {
        g <- graph.adjacency(edges, mode = "undirected", diag = FALSE)
    }
    plot(g,
       , layout       = layout.fruchterman.reingold(g)
       , edge.color   = "gray50"
       , vertex.color = "red"
       , vertex.size  = 3
       , vertex.label = NA
       , main         = main)
}




pdf("file.pdf", width = 8, height = 3)
par(mfrow = c(1,3), mar = c(2, 2, 2, 2))
data.plot(g, g.bic$index, sprintf("bic (%s sparsity)", round(g$sparsity[g.bic$index], 2)))
data.plot(g, g.ss$opt.index, sprintf("ss (%s sparsity)", round(g$sparsity[g.ss$opt.index], 2)))
data.plot(g, g.cv.min, sprintf("cv (%s sparsity)", round(g$sparsity[g.cv$minimizer], 2)))
dev.off()


sp <- get.sparsity(bma.edges)
pdf("file2.pdf", width = 8/3, height = 3)
par(mfrow = c(1,1), mar = c(2, 2, 2, 2))
data.plot(edges = bma.edges, main = sprintf("bma (%s sparsity)", round(sp, 2)))
dev.off()
