library(gRbase)
library(gRapHD)
library(mlbench)
library(glasso)

data(breastcancer)
b <- breastcancer[,1:1000] # High-dimensional: k=1000 >> n=249

c <- cor(b) # Empirical correlation matrix
d <- glasso(c, rho=.9) # LASSO with regularization parameter .9
d$loglik
g <- d$wi != 0 # Get edges
diag(g) <- F # Exclude self-edges 
rownames(g) <- colnames(g) <- names(breastcancer)[1:1000] # Get names of vertices
h <- as(g, "graphNEL")
p <- as(h, "gRapHD")
plot(p, numIter=1000)

system.time(d <- glasso(c, rho=.1))
system.time(d <- glasso(c, rho=.3))
system.time(d <- glasso(c, rho=.6))
system.time(d <- glasso(c, rho=.9))

