get.neighborhood <- function(edges, model.index, l) {
    model.edges <- edges$n.edges[model.index]
    ind <- setdiff(which(edges$n.edges <= model.edges + l &
                           edges$n.edges >= model.edges - l), model.index)
    if (length(ind) == 0) {
        #ind <- c(model.index + 1, model.index - 1)
        stop("No neighbors")
    }
        
    return(ind)
}

and.rule <- function(x) {
    for (i in 1:dim(x)[1]) 
        for (j in 1:dim(x)[1]) 
            if (x[i,j] != 1 | x[j,i] != 1) #(!all(x[i,j], x[j,i]))
                x[i,j] <- x[j,i] <- 0
    return(x)
}

get.edges.huge <- function(g) {
    n.edges <- simplify2array(lapply(g$path, sum)) / 2
    return(list(
        edges   = g$path,
        n.edges = n.edges,
        lambda  = g$lambda
        ))
}

ebic <- function(loglike, k, n) {
    -2*loglik + k*log(n) + 2*n*log(n)
}

bic <- function(loglik, k, n) {
    -n*loglik + k*log(n)
}


bmar <- function(x, l, niter, burnin) {
    n <- nrow(x); d <- ncol(x)
    THETA <- array(1, dim = c(d, d, niter))
    ## Initialize model space with regularization path
    #lambda.seq <- rev(exp(seq(log(0.05), log(0.2), length = 100)))
    g <- huge(x, nlambda = 1000, method = "glasso")
    lambda.seq <- g$lambda
            
    ## Obtain edges for all models in model space
    e <- get.edges.huge(g)
    stopifnot(all(e$n.edges == floor(g$df)))
    bic <- bic(g$loglik, e$n.edges, n)

    ## Initialize t=0 model with minimum sparsity
    model.t <- length(lambda.seq) #index

    accept <- 0
    for (t in 1:niter) {
        ## Propose new model on uniform distribution of neighbors
        ne <- get.neighborhood(e, model.t, l)
        model.new <- sample(ne, 1) #index

        ## Calculate MH acceptance ratio
        ## logB = log p(D|M1) / p(D|M2) approx= BIC1 - BIC2 / 2
        logr <- (bic[model.new] - bic[model.t]) / 2
        
        if (log(runif(1)) < logr) {
            model.t <- model.new
            accept <- accept + 1 # cannot deal with burnin
        }
        
        ## Compute theta hat (refit lasso)
        #THETA[, , t] <- glasso(s, rho = 0, zero = zero.index)$wi
        THETA[, , t] <- as.array(g$icov[[model.t]])
    }

    theta.hat <- apply(THETA, c(1, 2), mean)
    
    bmar <- list(
        theta.hat = theta.hat
      , l = l
      , niter = niter
      , burnin = burnin
      , accept = accept)
    class(bmar) <- "bmar"
    return(bmar)
}


