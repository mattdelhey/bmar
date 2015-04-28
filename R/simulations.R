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
        d <- huge(x.trn, lambda = lambda, method = "glasso")
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
    for (i in 1:length(lambda.seq)) {
        message(sprintf("[%i of %i] \t lambda = %f", i, length(lambda.seq), lambda.seq[i]))
        l <- g.folds(K=K, x=x, s=s, lambda=lambda.seq[i])
        n.edges[i] <- l$n.edges
        m[i] <- l$mean
        sd[i] <- l$sd
    }
    minimizer <- which.min(m)
    return(list(minimizer = minimizer, n.edges = n.edges, m = m, sd = sd))
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
    c(sse       = g.sse(theta, theta.hat)
    , kl        = g.kl(theta, theta.hat, sigma, d)
    , precision = g.precision(edges, edges.hat)
    , recall    = g.recall(edges, edges.hat, d)
    , f1        = g.f1(edges, edges.hat, d))
}

bic <- function(loglik, k, n) {
    -n*loglik + k*log(n)
}

select.bic <- function(g, n, d) {
    n.edges <- simplify2array(lapply(g$path, sum))
    scores <- bic(g$loglik, n = n, k = d*d - n.edges)
    index <- which.min(scores)
    return(list(scores = scores, index = index))
}

sim.summarize <- function(sim.objects, nlambda, K, stars.thresh, stars.subsample.ratio, rep.num) {
    simplify2array(lapply(sim.objects, function(sim) {        
        theta <- sim$omega
        sigma <- sim$sigma
        edges <- sim$theta
        x <- sim$data; d <- ncol(x); n <- nrow(x)
        ## Model estimation via GLASSO
        message(sprintf("[simulation] \t Simulation for graph type {%s}", sim$graph.type))
        g <- huge(x, nlambda = nlambda, method = "glasso")
        lambda.seq <- g$lambda
        # Model selection methods
        message(sprintf("... [method] \t BIC..."))
        g.bic <- select.bic(g = g, n = n, d = d)
        message(sprintf("... [method] \t SS..."))
        g.ss <- huge.select(g, criterion = "stars",
                            stars.thresh = stars.thresh,
                            stars.subsample.ratio = stars.subsample.ratio,
                            rep.num = rep.num)
        message(sprintf("... [method] \t CV..."))
        g.cv <- g.cvf(K = K, x = x, lambda.seq = lambda.seq)
        ## Sparsity
        sparsity <- g$sparsity[c(g.bic$index, g.ss$opt.index, g.cv$minimizer)]
        message(sprintf("[sparsity] %s \t :: %s \n", c("BIC", "SS", "CV"), round(sparsity, 2)))
        ## Error
        error.bic <- error.functions(g$icov[[g.bic$index]], g$path[[g.bic$index]], theta, sigma, edges, d)
        error.ss  <- error.functions(g$icov[[g.ss$opt.index]], g.ss$refit,  theta, sigma, edges, d)
        error.cv  <- error.functions(g$icov[[g.cv$minimizer]], g$path[[g.cv$minimizer]], theta, sigma, edges, d)                          
        return(rbind(
            bic = c(error.bic, sparsity = sparsity[1]),
            ss  = c(error.ss,  sparsity = sparsity[2]),
            cv  = c(error.cv,  sparsity = sparsity[3])))        
    }))
}
    
