UniqInd <- function(z)
{
    z <- as.integer(z)
    z.sorted <- sort.int(z, index.return = TRUE)
    sortZ <- z.sorted$x
    indSortZ <- z.sorted$ix
    numelZ <- length(z)
    
    if ((is.numeric(sortZ) && (numelZ>1)))
    {
        dSortZ <- diff(sortZ)
        if ((is.nan(dSortZ[1])) || (is.nan(dSortZ[numelZ-1])))
            groupsSortZ <- sortZ[1:numelZ-1] != sortZ[2:numelZ]
        else groupsSortZ <- dSortZ != 0
    }else groupsSortZ <- sortZ[1:numelZ-1] != sortZ[2:numelZ]
    
    if (numelZ!=0)
    {
        groupsSortZ <- c(TRUE, groupsSortZ)
    }else groupsSortZ <- matrix(0,0,0)
    
    ## Find indices
    indC <- cumsum(as.numeric(groupsSortZ))
    indC[indSortZ] <- indC
    
    return(indC)
}