#' Gamma-Negative Binomial process modeling of sequencing count data
#'
#' This function fits GNBP model to counts in different conditions and returns the
#' posterior samples of probability and dispersion parameters than can be used to
#' assess the differential expression significance.
#'
#' @param y a matrix of counts, rows are corresponding to genes and columns are
#' corresponding to samples.
#' @param groups a vector of factors indicating group memberships of samples
#' @param randtry state of set.seed
#' @param Burnin Number of burn-in iterations in MCMC
#' @param Collections Number of collected posterior samples after burn-in
#' @return posterior samples of probability and dispersion parameters of Negative Binomial
#' @examples
#' y <- rnbinom(n = 100, size = 5, prob = 0.5)
#' y <- matrix(data = y, nrow = 10, ncol = 10)
#' groups <- factor(rep(c(1,2), each = 5))
#' res <- GNBP(y, groups, randtry = 1, Burnin = 1000L, Collections = 1000L)
#' r <- res$r_k
#' gene.kl <- KLsym(r[,1:1000], r[,1001:2000])
#' @export
GNBP <- function(y, groups, randtry=1, Burnin=1000L, Collections=1000L)
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
    r <- matrix(1, nrow=ngenes, ncol=ngroups)
    r_star <- rep(0, ngroups)
    p <- rep(0.5, nsamples)
    gamma0 <- rep(1, ngroups)
    c <- rep(1, ngroups)

    # Hyperparameters
    a0 <- b0 <- e0 <- f0 <- c0 <- d0 <- 1e-3
    realmin <- 2.2251e-308

    # Loop indices
    counter <- 1L
    iterMax <- Burnin+Collections

    post.r <- matrix(0, nrow = ngenes, ncol = ngroups*Collections/A)
    post.p <- matrix(0, nrow = nsamples, ncol = Collections/A)
   # r.tilde <- matrix(0, nrow = length(idx.test), ncol = ngroups*Collections/A)
    for (iter in 1:iterMax)
    {
   #     cat(iter, '\n')
        for (g in glevels)
        {
            subg <- which(groups==g)
            Kj <- length(idx.nz[[g]])
            idx <- idx.nz[[g]]
            q.dot <- sum(log(pmax(1 - p[subg], realmin)))
            p.prime <- -q.dot/(c[g] - q.dot)
            gamma0[g] <- rgamma(1, e0 + Kj, rate = f0 - log(max(1-p.prime, realmin)))
            L <- rep(0, Kj)
            for (j in subg)
            {
                L <- L + CRT_vector(y[idx, j],r[idx, g])
            }
            r[idx, g] <- rgamma(Kj, L, rate = c[g] - q.dot)
            r_star[g] <- rgamma(1, gamma0[g], rate = c[g] - q.dot)
            mass <- sum(r[idx, g]) + r_star[g]
            c[g] <- rgamma(1, c0 + gamma0[g], rate = d0 + mass)
            p[subg] <- rbeta(length(subg), a0 + lib.size[subg], b0 + mass)

            if (iter>Burnin)
            {
                post.r[idx, (g-1)*Collections + (iter-Burnin)] <- r[idx, g]
             #   r.tilde[, (g-1)*Collections/A + (iter-Burnin)/A] <- post.r[idx.test, (g-1)*Collections/A + (iter-Burnin)/A] +
             #                                                   rgamma(length(idx.test), gamma0[g]/ngenes, rate = c[g] - q.dot)
            }
        }
    if (iter>Burnin)  post.p[,(iter-Burnin)] <- p
    }
    return(list(r_k = post.r, gamma0 = gamma0, p_j = post.p))
}
