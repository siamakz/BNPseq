#' Scaled Negative Binomial process modeling of sequencing count data
#'
#' This function fits scaled NBP model to counts in different conditions and returns the
#' posterior samples of normalized rate parameters than can be used to
#' assess the differential expression significance.
#'
#' @param y a matrix of counts, rows are corresponding to genes and columns are
#' corresponding to samples.
#' @param groups a vector of factors indicating group memberships of samples
#' @param randtry state of set.seed
#' @param Burnin Number of burn-in iterations in MCMC
#' @param Collections Number of collected posterior samples after burn-in
#' @return posterior samples of rate parameters of Negative Binomial
#' @examples
#' y <- rnbinom(n = 100, size = 5, prob = 0.5)
#' y <- matrix(data = y, nrow = 10, ncol = 10)
#' groups <- factor(rep(c(1,2), each = 5))
#' res <- NBP(y, groups, randtry = 1, Burnin = 1000L, Collections = 1000L)
#' r <- res$r
#' gene.kl <- KLsym(r[,1:1000], r[,1001:2000])
#' @export
NBPscaled <- function(y, groups, randtry=1, Burnin=500L, Collections=1000L)
{
    set.seed(randtry)

    y <- as.matrix(y)
    groups <- as.factor(groups)
    ngenes <- nrow(y)
    nsamples <- ncol(y)
    lib.size <- colSums(y)
    glevels <- as.integer(levels(groups))
    ngroups <- length(glevels)
    idx.nz <- vector("list", ngroups)
    for (g in glevels)
    {
        subg <- which(groups==g)
        idx.nz[[g]] <- which(rowSums(y[, subg])!=0)
    }

    # Initialize parameters
    r <- matrix(0, nrow=ngenes, ncol=ngroups)
    r_star <- rep(0, ngroups)
    gamma0 <- rep(1, ngroups)
    c <- rep(1, ngroups)
    q <- rep(1, nsamples)

    # Hyperparameters
    a0 <- b0 <- e0 <- f0 <- c0 <- d0 <- 1e-3
    realmin <- 2.2251e-308

    # Loop indices
    counter <- 1L
    iterMax <- Burnin+Collections

    post.r <- matrix(0, nrow = ngenes, ncol = ngroups*Collections)
    for (iter in 1:iterMax)
    {
        cat(iter, '\n')
        for (g in glevels)
        {
            subg <- which(groups==g)
            J <- length(subg)
            Kj <- length(idx.nz[[g]])
            q.dot <- sum(q[subg])
            idx <- idx.nz[[g]]
            gamma0[g] <- rgamma(1, e0 + Kj, rate = f0 - log(max(c[g]/(c[g]+q.dot), realmin)))
            r[idx, g] <- rgamma(Kj, rowSums(y[idx,subg]), rate = c[g]+q.dot)
            r_star[g] <- rgamma(1, gamma0[g], rate = c[g]+q.dot)
            mass <- sum(r[idx, g]) + r_star[g]
            c[g] <- rgamma(1, c0 + gamma0[g], rate = d0 + mass)
            q[subg] <- rgamma(J, a0+lib.size[subg], rate = b0 + mass)

            if (iter>Burnin)
            {
                post.r[idx, (g-1)*Collections + iter-Burnin] <- r[idx, g]
            }
        }
    }
    return(list(r = post.r, gamma0 = gamma0))
}
