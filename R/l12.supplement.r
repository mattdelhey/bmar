library(glasso)
library(gRapHD)
library(cvTools)
library(gRain)
library(spcov)
library(igraph)
library(huge)

cv <- function(K, x, s, lambda)
{
  f <- vector(length=K)
  for (i in 1:K)
    {
    x.training <- x[s$which!=i,]
    d <- glasso(cov(x.training), rho=lambda, penalize.diagonal=FALSE)
    g <- d$wi != 0
    diag(g) <- FALSE
    number.edges <- sum(g==TRUE)
    x.testing <- x[s$which==i,]
    c.testing <- cov(x.testing)
    f[i] <- - log(det(d$wi)) + sum(diag(c.testing %*% d$wi)) 
    }
  l <- list()
  l$number.edges <- number.edges
  l$m <- mean(f)
  l$sd <- sd(f)
  l
}

# Generating uncorrelated gaussian data or sparse graph with cliques of size 2
library(MASS)
t <- 100
n <- 100
#Sigma <- diag(n)
Sigma <-  GenerateCliquesCovariance(ncliques=50, cliquesize=2, 1)$Sigma
x <- mvrnorm(n = t, rep(0, n), Sigma)

# Plot true graph

diag(Sigma) <- FALSE
g1 <- graph.adjacency(Sigma!=0)
layout(1)
plot(g1, vertex.size = 3,vertex.label = NA, edge.arrow.size=.2)

# Cross-validation

K <- 5
s <- cvFolds(n=nrow(x), K=K, type="random")
number.candidates <- 20 
number.edges <- vector(length=number.candidates)
m <- vector(length=number.candidates)
sd <- vector(length=number.candidates)
lambda <-  2 ^ (c(1:number.candidates) / (2*number.candidates)) - 1
for (i in 1:number.candidates) 
{
  print(lambda[i])
  l <- cv(K=K, x=x, s=s, lambda=lambda[i])
  number.edges[i] <- l$number.edges
  m[i] <- l$m
  sd[i] <- l$sd
}

# Output

plot(lambda, m, xlab=expression(lambda), ylab=expression(f), main="", type="b", col="blue")

# cv minimizer

minimizer <- which.min(m)

# Model selection using cv

c <- cov(x)
d <- glasso(c, rho=lambda[minimizer], penalize.diagonal = FALSE)
g <- d$wi != 0
diag(g) <- FALSE
g1 <- graph.adjacency(g!=0)
layout(1)
plot(g1, vertex.size = 3, vertex.label = NA, edge.arrow.size = .2)

# Model selection using stars

out.glasso <- huge(x, method = "glasso")
out.select <- huge.select(out.glasso, criterion = "stars", stars.thresh = 0.05, rep.num = 10)
plot(out.select)
