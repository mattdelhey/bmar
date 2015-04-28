#library(glasso)

## replace.diag <- function(x, replacement) {
##     K <- dim(x)[3]
##     for (k in 1:K) {
##         diag(x[, , k]) <- replacement
##     }
##     x
## }
## to.zero.matrix <- function(x) {
##     ind <- which(!x, arr.ind = TRUE)
##     matrix(ind[1:(nrow(ind)/2), ], ncol = 2)
## }
## get.edges2 <- function(path, diag = FALSE) {    
##     edges <- simplify2array(lapply(1:dim(path$wi)[3], function(k) {
##         e <- path$wi[, , k] != 0
##         diag(e) <- FALSE
##         e <- and.rule(e)
##     }))
##     n.edges <- apply(edges, 3, function(i) length(which(i)) / 2)
##     return(list(
##         edges   = edges
##       , n.edges = n.edges
##       , lambda  = path$rholist))    
## }


get.neighborhood <- function(edges, model.index, l) {
    model.edges <- edges$n.edges[model.index]
    ind <- setdiff(which(edges$n.edges <= model.edges + l &
                           edges$n.edges >= model.edges - l), model.index)
    if (length(ind) == 0)
        ind <- c(model.index + 1, model.index - 1)
    return(ind)
}

and.rule <- function(x) {
    for (i in 1:dim(x)[1]) 
        for (j in 1:dim(x)[1]) 
            if (x[i,j] != 1 | x[j,i] != 1) #(!all(x[i,j], x[j,i]))
                x[i,j] <- x[j,i] <- 0
    return(x)
}

get.edges <- function(path) {   
    edges <- lapply(path$path, function(m) {
        m <- as.matrix(m)
        diag(m) <- 0
        m <- and.rule(m)
    })

    #raw.edges <- simplify2array(lapply(path$path, function(m) length(slot(m, "i"))))
    #all(simplify2array(lapply(edges, isSymmetric)))

    n.edges <- simplify2array(lapply(edges, function(i) length(which(i == 1)))) / 2   

    return(list(
        edges   = simplify2array(edges)
      , n.edges = n.edges
      , lambda  = path$lambda))
}

ebic <- function(loglike, k, n) {
    -2*loglik + k*log(n) + 2*n*log(n)
}

l <- 250
accept <- 0
bmar <- function(x, l, niter, burnin) {
    n <- nrow(x); d <- ncol(x)
    theta.hat <- array(1, dim = c(d, d, niter))
    ## Initialize model space with regularization path
    lambda.seq <- rev(exp(seq(log(0.05), log(0.2), length = 100)))
    path <- huge(x, nlambda = 1000, method = "glasso")
    lambda.seq <- path$lambda
            
    ## Obtain edges for all models in model space
    e <- get.edges(path)
    stopifnot(all(e$n.edges == floor(path$df)))
    bic <- get.bic(path$loglik, e$n.edges, n)

    ## Initialize t=0 model with minimum sparsity
    model.t.index <- length(lambda.seq)
    
    for (t in 1:niter) {
        ## Propose new model on uniform distribution of neighbors
        ne <- get.neighborhood(e, model.t.index, l)
        model.new.index <- sample(ne, 1)

        ## Calculate MH acceptance ratio
        ## logB = log p(D|M1) / p(D|M2) approx= BIC1 - BIC2 / 2
        logr <- bic[model.new.index] - bic[model.t.index] / 2
        
        if (log(runif(1)) < logr) {
            model.t.index <- model.new.index
            accept <- accept + 1 # cannot deal with burnin
        }
        
        ## Compute theta hat (refit lasso)
        zero.index <- to.zero.matrix(e$edges[, , model.new.index])
        theta.hat[, , t] <- glasso(s, rho = 0, zero = zero.index)$wi

        ## Generate zero constraints for new model
        #zero.index <- to.zero.matrix(e$edges[, , model.new.index])
        #stopifnot(nrow(zero.index) == d^2 - e$n.edges[model.new.index])
        
        #model.new <- glasso(s, rho = 0, zero = zero.index, start = "warm",
        #                    w.init  = path$w[, , model.new.index],
        #                    wi.init = path$wi[, , model.new.index])

    }

    theta.hat <- apply(theta.hat, c(1, 2), mean)
    
    return(list(
        theta.hat = theta.hat
      , l = l
      , niter = niter
      , burnin = burnin))
}

bmar(x, l=3, niter=10, burnin=0)

