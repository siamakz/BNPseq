logBetarnd <- function(n, gamma0, cc)
{
    tol <- 1e-7
    x0 <- 1
    xinc <- 2
    m <- 11
    L <- 1
    A <- 19
    nburn <- 38
    maxiter <- 500
    
    nterms <- nburn + m*L
    seqbtL <- seq(nburn, nterms, L)
    y <- pi * (1i) * (1:nterms) / L
    expy <- exp(y)
    A2L <- 0.5 * A / L
    expxt <- exp(A2L) / L
    
    coef <- exp( lgamma(m+1)-lgamma(1:(m+1))-lgamma(seq(m+1, 1,by=-1)) -m*log(2) )
    u <- sort(runif(n))
    xrand <- u
    
    t <- x0/xinc
    cdf <- 0
    kount0 <- 0
    set1st <- FALSE
    while (kount0 < maxiter && cdf < u[n])
    {
        t <- xinc * t
        kount0 <- kount0 + 1
        x <- A2L / t
        z <- x + y/t
        ltx <- ltpdf(x, gamma0,cc)
        ltzexpy <- ltpdf(z, gamma0,cc)*expy
        par.sum <- 0.5*Re(ltx) + cumsum( Re(ltzexpy) )
        par.sum2 <- 0.5*Re(ltx/x) + cumsum( Re(ltzexpy/z) )
        pdf <- expxt*sum(coef* par.sum[seqbtL]) / t
        cdf <- expxt*sum(coef* par.sum2[seqbtL]) / t
        if (!set1st && cdf > u[1])
        {
            cdf1 <- cdf
            pdf1 <- pdf
            t1 <- t
            set1st <- TRUE
        }
    }
    if (kount0 >= maxiter)
    {
        stop('Cannot locate upper quantile')
    }
    
    upplim <- t
    
    lower <- 0
    t <- t1
    cdf <- cdf1
    pdf <- pdf1
    kount <- rep(0,n)
    
    maxiter <- 1000
    
    for (j in 1:n)
    {
        upper <- upplim
        kount[j] <- 0
        while (kount[j] < maxiter && abs(u[j]-cdf) > tol)
        {
            kount[j] <- kount[j] + 1
            t <- t - (cdf-u[j])/pdf
            if (t < lower || t > upper)
            {
                t <- 0.5 * (lower + upper)
            }
            
            x <- A2L / t
            z <- x + y/t
            ltx <- ltpdf(x, gamma0,cc)
            ltzexpy <- ltpdf(z, gamma0,cc)* expy
            par.sum <- 0.5*Re(ltx) + cumsum( Re(ltzexpy) )
            par.sum2 <- 0.5*Re(ltx/x) + cumsum( Re(ltzexpy/z) )
            pdf <- expxt * sum(coef * par.sum[seqbtL]) / t
            cdf <- expxt * sum(coef * par.sum2[seqbtL]) / t
            
            if (cdf <= u[j]) lower <- t else upper <- t
        }
        
        if (kount[j] >= maxiter) cat('Desired accuracy not achieved for F(x)<-u')
        
        xrand[j] <- t
        lower <- t
    }
    
    if (n > 1)
        rsample <- xrand[sample(n)]
    else
        rsample <- xrand
    
    rsample
}