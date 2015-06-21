##' Find the significant eigenvalues of a matrix.
##'
##' @title Tracy-Wisdom test
##' @param eigenvalues a numeric vector whose elements are the
##' eigenvalues of a matrix. The values should be sorted in the
##' descending order.
##' @param eigenL the number of the eigenvalues.
##' @param criticalpoint a numeric value corresponding to the
##' significance level. If the significance level is 0.05, 0.01,
##' 0.005, or 0.001, the critical point should be set to be 0.9793,
##' 2.0234, 2.4224, or 3.2724, accordingly. The default is 2.0234.
##' @return A list of \code{SigntEigenL} and
##' \code{TW.stat}. \code{SigntEigenL} is the number of the
##' significant eigenvalues and \code{TW.stat} is a vector of the
##' Tracy-Wisdom statistics.
##' @author Lin Wang, Wei Zhang, and Qizhai Li.
##' @references N Patterson, AL Price, and D Reich. Population
##' structure and eigenanalysis. \emph{PloS Genetics}. 2006; 2(12):
##' 2074-2093.
##' @references CA Tracy and H Widom. Level-spacing distributions and
##' the Airy kernel. \emph{Communications in Mathematical
##' Physics}. 1994; 159(1): 151-174.
##' @examples
##' tw(eigenvalues = c(5, 3, 1, 0), eigenL = 4, criticalpoint = 2.0234)
##' @export
tw <- function(eigenvalues, eigenL, criticalpoint=2.0234)
{
    dex <- which(eigenvalues <= 1e-8)
    eigenvalues[dex] <- 1e-8

    L1 <- rev(cumsum(rev(eigenvalues)))
    L2 <- rev(cumsum(rev(eigenvalues^2)))
    N <- eigenL:1
    S2 <- N^2*L2/(L1^2)
    v <- N*(N+2)/(S2-N) # Effective number of markers

    L <- N*eigenvalues/L1

    v.st <- sqrt(v-1)
    N.st <- sqrt(N)

    mu  <- (v.st+N.st)^2/v
    sig <- (v.st+N.st)/v * (1/v.st+1/N.st)^(1/3)

    twstat <-(L-mu)/sig

    #sink(output)
    #cat("TWstat = ", twstat, '\n')
    #sink()

    d <- which(twstat < criticalpoint)[1]

    if (length(d)==0)
    {
        d <- -100
    }else
    {
        d <- d-1
    }

    res <- list(SigntEigenL=d, TW.stat=twstat)

    return(res)
}
