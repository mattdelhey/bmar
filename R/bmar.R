get.neighborhood <- function(edges, model.index, l) {
    model.edges <- edges$n.edges[model.index]
    max.index <- length(edges$n.edges)
    ind <- setdiff(which(edges$n.edges <= model.edges + l &
                           edges$n.edges >= model.edges - l), model.index)
    if (length(ind) == 0) {
        #warning("No neighbors")
        ind <- c(min(max.index, model.index + 1), max(1, model.index - 1))
    }
        
    return(ind)
}

and.rule <- function(x) {
    for (i in 1:dim(x)[1]) {
        for (j in 1:dim(x)[2]) {
            if (x[i,j] != 0 & x[j,i] != 0) {
                #(!all(x[i,j], x[j,i]))
                x[i,j] <- x[j,i] <- 1
            } else {
                x[i,j] <- 0
            }
        }
    }                
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

sparse.prior <- function(edges, d) {
    #m.l1 <- sum(abs(m))
    m.l1 <- sum(edges)
    (  m.l1 / ( exp(1)*d*(d-1) )  )^(m.l1)
}


bmar <- function(x = NULL, g = NULL, l, nlambda = 1000, niter, burnin) {
    ## Initialize model space with regularization path
    if (!is.null(x)) {
        n <- nrow(x); d <- ncol(x)
        #lambda.seq <- rev(exp(seq(log(0.05), log(0.2), length = 100)))
        g <- huge(x, nlambda = nlambda, method = "glasso")        
    }
    if (!is.null(g)) {
        n <- nrow(g$data); d <- ncol(g$data)
    }

    lambda.seq <- g$lambda    
    THETA <- array(1, dim = c(d, d, niter))        
            
    ## Obtain edges for all models in model space
    e <- get.edges.huge(g)
    #stopifnot(all(e$n.edges == floor(g$df))) #these will differ
    bic <- bic(g$loglik, e$n.edges, n)

    ## Initialize t=0 model with minimum sparsity
    model.t <- length(lambda.seq) #index

    accept <- 0
    for (t in 1:(niter+burnin)) {
        ## Propose new model on uniform distribution of neighbors
        ne <- get.neighborhood(e, model.t, l)
        model.new <- sample(ne, 1) #index

        ## Calculate MH acceptance ratio
        ## logB = log p(D|M1) / p(D|M2) approx= BIC1 - BIC2 / 2
        # BIC approximation
        #logr <- (bic[model.new] - bic[model.t]) / 2
        # Flat prior
        logr <- g$loglik[model.new] - g$loglik[model.t]
        # Sparsity prior
        #logr <- ( g$loglik[model.new] + sparse.prior(g$path[[model.new]], d) ) -
          #( g$loglik[model.t] + sparse.prior(g$path[[model.t]], d) )
        
        if (log(runif(1)) < logr) {
            model.t <- model.new
            accept <- accept + 1 # cannot deal with burnin
        }
        
        ## Compute theta hat (refit lasso)
        if (t > burnin)
            THETA[, , t-burnin] <- as.array(g$icov[[model.t]])
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


