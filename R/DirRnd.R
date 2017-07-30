DirRnd <- function(alpha)
## Dirichlet random numbers generator
{
    y <- rgamma(length(alpha), alpha)
    return(y/sum(y))
}