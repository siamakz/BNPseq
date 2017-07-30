ltpdf <- function(s, gamma0, cc)
{
    x <- exp(-gamma0*(psi_complex(cc+s)-psi_complex(cc)))
    x
}