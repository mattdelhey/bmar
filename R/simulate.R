#' get.lambda.seq
#' Create an appropriately scaled sequence of candidate lambda values
#' @param S correlation of data
get.lambda.seq <- function(nlambda, S, d, lambda.min.ratio = 0.1) {
    lambda.max <- max(max(S-diag(d)),-min(S-diag(d)))
    lambda.min <- lambda.min.ratio*lambda.max
    lambda <- exp(seq(log(lambda.max), log(lambda.min), length = nlambda))
}

redu <- function(x, l) {
    ind <- floor(seq(1, length(x), length.out = l))
    x[ind]
}

get.sparsity <- function(edges) {
    n <- prod(dim(edges))
    n.nz <- length(which(edges == 1))
    n.nz / n
}

ose.rule <- function(m, sd) {
    min <- which.min(m)
    sd.min <- sd[min]
    which.max(m <= m[min] + sd.min)
}
    
g.folds <- function(K, x, s, lambda) {
    stopifnot(length(lambda) == 1)
    f <- vector(length = K)
    for (k in 1:K) {
        x.trn <- x[s$which != k, ]
        d <- huge(x.trn, lambda = lambda, method = "glasso", verbose = FALSE)
        theta <- as.matrix(d$icov[[1]])
        edges <- as.matrix(d$path[[1]])
        n.edges <- sum(edges)
        x.tst <- x[s$which == k, ]
        f[k] <- -log(det(theta)) + tr(cov(x.tst) %*% theta)
    }
    return(list(
        n.edges = n.edges
      , mean    = mean(f)
      , sd      = sd(f)))
}

g.cvf <- function(K, x, lambda.seq) {
    s <- cvTools::cvFolds(n=nrow(x), K=K, type="random")
    m <- sd <- n.edges <- vector(length = length(lambda.seq))

    library(parallel)
    cl <- makeCluster(4)
    clusterEvalQ(cl, {library(huge); library(bmar)})
    
    r <- t(simplify2array(parLapply(cl, 1:length(lambda.seq), function(i) {
        mes <- sprintf("Conducting cross-validation....in progress: [%i of %i] \t lambda = %f", i, length(lambda.seq), lambda.seq[i])
        cat(mes, "\r")
        flush.console()        
        l <- g.folds(K=K, x=x, s=s, lambda=lambda.seq[i])        
        c(n.edges = l$n.edges
        , m = l$mean
        , sd = l$sd)
    })))
    
    stopCluster(cl)

    minimizer <- which.min(r[, "m"])
    return(list(minimizer = minimizer, n.edges = r[, "n.edges"],
                m = r[, "m"], sd = r[, "sd"]))
}

sim.plot <- function(sim.objects, f.out) {
    pdf(f.out, width = 6, height = 6)
    par(mfrow = c(2,2), mar = c(2, 2, 2, 2))
    lapply(sim.objects, function(sim) {
        g <- graph.adjacency(as.matrix(sim$theta), mode = "undirected", diag = FALSE)
        plot(g,
           , layout       = layout.fruchterman.reingold(g)
           , edge.color   = "gray50"
           , vertex.color = "red"
           , vertex.size  = 3
           , vertex.label = NA
           , main         = sim$graph.type)
    })
    dev.off()
}

g.sse <- function(theta, theta.hat) {
    norm(theta.hat - theta, type = "F")^2
}

tr <- function(x) {
    sum(diag(x))
}

g.kl <- function(theta, theta.hat, sigma, d) {
    -log(det(theta.hat)) + tr(theta.hat %*% sigma) - (-log(det(theta)) + d)
}

g.precision <- function(edges, edges.hat) {
    tp.all <- (edges != 0) * (edges.hat !=0)
    diag(tp.all) <- 0
    sum(tp.all != 0) / sum(edges != 0)
}

g.recall <- function(edges, edges.hat, d) {
    fp.all <- (edges == 0) * (edges.hat != 0)
    diag(fp.all) <- 0
    sum(fp.all != 0) / (d*(d-1) - sum(edges.hat != 0))
}

g.f1 <- function(edges, edges.hat, d) {
    precision <- g.precision(edges, edges.hat)
    recall <- g.recall(edges, edges.hat, d)
    2*(precision * recall) / (precision + recall)
}

error.functions <- function(theta.hat, edges.hat, theta, sigma, edges, d) {
    theta <- as.matrix(theta)
    theta.hat <- as.matrix(theta.hat)
    c(sse       = g.sse(theta, theta.hat)
    , kl        = g.kl(theta, theta.hat, sigma, d)
    , precision = g.precision(edges, edges.hat)
    , recall    = g.recall(edges, edges.hat, d)
    , f1        = g.f1(edges, edges.hat, d))
}

select.bic <- function(g, n, d) {
    n.edges <- simplify2array(lapply(g$path, sum))
    scores <- bic(g$loglik, n = n, k = d*d - n.edges)
    index <- which.min(scores)
    return(list(scores = scores, index = index))
}

sim.summarize <- function(sim.objects, nlambda, K, stars.thresh,
                          stars.subsample.ratio, rep.num,
                          l, niter, burnin) {
    simplify2array(lapply(sim.objects, function(sim) {        
        theta <- as.matrix(sim$omega)
        sigma <- sim$sigma
        edges <- sim$theta
        x <- sim$data; d <- ncol(x); n <- nrow(x)
        
        ## Model estimation via GLASSO
        message(sprintf("[simulation] \t Simulation for graph type {%s}", sim$graph.type))
        g <- huge(x, nlambda = nlambda, method = "glasso")
        lambda.seq <- g$lambda
        lambda.seq.cv <- redu(g$lambda, length(g$lambda)/10)
        
        
        ## Model selection methods
        
        # BMA
        message(sprintf("... [method] \t BMA..."))
        t.bma <- system.time(g.bma <- bmar(g = g, l = l, niter = niter, burnin = burnin))[3]
        bma.edges <- and.rule(g.bma$theta.hat)

        # BIC
        message(sprintf("... [method] \t BIC..."))
        t.bic <- system.time(g.bic <- select.bic(g = g, n = n, d = d))[3]

        # SS
        message(sprintf("... [method] \t SS..."))
        t.ss <- system.time(g.ss <- huge.select(g, criterion = "stars",
                            stars.thresh = stars.thresh,
                            stars.subsample.ratio = stars.subsample.ratio,
                            rep.num = rep.num))[3]
        # CV
        message(sprintf("... [method] \t CV..."))
        t.cv <- system.time(g.cv <- g.cvf(K = K, x = x, lambda.seq = lambda.seq.cv))[3]
        
        ## Sparsity
        sparsity <- c(g$sparsity[c(g.bic$index, g.ss$opt.index, g.cv$minimizer)], get.sparsity(bma.edges))
        message(sprintf("[sparsity] %s \t :: %s \n", c("BIC", "SS", "CV", "BMA"), round(sparsity, 2)))
        
        ## Error
        error.bic <- error.functions(g$icov[[g.bic$index]], g$path[[g.bic$index]], theta, sigma, edges, d)
        error.ss  <- error.functions(g.ss$opt.icov, g.ss$refit,  theta, sigma, edges, d)
        error.cv  <- error.functions(g$icov[[g.cv$minimizer]], g$path[[g.cv$minimizer]], theta, sigma, edges, d)
        error.bma <- error.functions(g.bma$theta.hat, bma.edges, theta, sigma, edges, d)
        return(rbind(
            bic = c(error.bic, sparsity = sparsity[1], t.bic),
            ss  = c(error.ss,  sparsity = sparsity[2], t.ss),
            cv  = c(error.cv,  sparsity = sparsity[3], t.cv),
            bma = c(error.bma, sparsity = get.sparsity(bma.edges), t.bma)))
    }))
}
    
