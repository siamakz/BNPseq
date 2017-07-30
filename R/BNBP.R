#' Beta-Negative Binomial process modeling of sequencing count data
#'
#' This function fits BNBP model to counts in different conditions and returns the
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
#' res <- BNBP(y, groups, randtry = 1, Burnin = 1000L, Collections = 1000L)
#' p <- res$p_k
#' p <- p/(1-p)
#' gene.kl <- KLsym(p[,1:1000], p[,1001:2000])
#' @export
BNBP <- function(y, groups, randtry = 1, Burnin = 1000L, Collections = 1000L)

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
    p <- matrix(1, nrow=ngenes, ncol=ngroups)
    p_star <- rep(0, ngroups)
    r <- rep(1, nsamples)
    gamma0 <- rep(1, ngroups)
    c <- rep(1, ngroups)

    # Hyperparameters
    a0 <- b0 <- e0 <- f0 <- 1e-3
    c0 <- d0 <- 1
    realmin <- 2.2251e-308

    # Loop indices
    counter <- 0L
    iterMax <- Burnin+Collections

    post.p <- matrix(0, nrow = ngenes, ncol = ngroups*Collections)
    post.r <- matrix(0, nrow = nsamples, ncol = Collections)
    for (iter in 1:iterMax)
    {
        cat(iter, '\n')
        for (g in glevels)
        {
            subg <- which(groups==g)
            Kj <- length(idx.nz[[g]])
            idx <- idx.nz[[g]]
            r.dot <- sum(r[subg])
            gamma0[g] <- rgamma(1, e0 + Kj, rate = f0 + digamma(c[g]+r.dot) - digamma(c[g]))
            p[idx, g] <- rbeta(Kj, rowSums(y[idx, subg]), c[g]+r.dot)
            p_star[g] <- logBetarnd(1, gamma0[g], c[g]+r.dot)

            #L <- CRT_matrix(y[, subg],r[subg])
            L <- rep(0, length(subg))
            for (k in idx)
            {
                L <- L+CRT_vector(y[k, subg], r[subg])
            }
            sumlogpk <- sum(log(pmax(1-p[idx,g], realmin)))
            r[subg] <- rgamma(length(subg), a0+L, rate = b0 + p_star - sumlogpk)
            # MH for c
            sumlogpk <- -sum(p[idx, g])
            r.dot <- sum(r[subg])
            cnew <- rgamma(1, c0+gamma0[g], 1/(d0 - sumlogpk + p_star[g]))

            temp <- (-gamma0[g]*(digamma(r.dot+cnew)-digamma(cnew)) +
                    Kj*lgamma(cnew+r.dot) - sum(lgamma(cnew+r.dot+rowSums(y[, subg])))  + 0*Kj*log(cnew) ) +
                (c0-1)*log(cnew)-d0*cnew
            temp <- temp - (-gamma0[g]*(digamma(r.dot+c[g])-digamma(c[g])) +
                          Kj*lgamma(c[g]+r.dot) - sum(lgamma(c[g]+r.dot+rowSums(y[, subg])))  + 0*Kj*log(c[g]))
            temp <- temp -( (c0-1)*log(c[g])-d0*c[g] )
            temp <- temp + (c0+gamma0[g]-1)*log(c[g]) - (d0 -sumlogpk + 1*p_star[g])*c[g]
            temp <- temp - ((c0+gamma0[g]-1)*log(cnew) - (d0 -sumlogpk + 1*p_star[g])*cnew)

            if (exp(temp)>runif(1))
            {
                c[g] <- cnew
                counter <- counter + 1
            }

            if (iter>Burnin)
            {
                post.p[idx, (g-1)*Collections + iter-Burnin] <- p[idx, g]
            }
        }
    if (iter>Burnin)
    {
        post.r[,iter-Burnin] <- r
    }
    }
    return(list(p_k = post.p, r_j = post.r))
}
