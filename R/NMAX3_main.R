##' Test for the association between a biallelic SNP and a
##' quantitative trait using the maximum value of the three
##' nonparametric trend tests derived for the recessive, additive, and
##' dominant models. It is a robust procedure against the genetic
##' models.
##'
##' Under the null hypothesis of no association, the vector of the
##' three nonparametric tests under the recessive, additive, and
##' dominant models asymptotically follows a three-dimensional normal
##' distribution. The p-value can be calculated using the function
##' \link{pmvnorm} in the R package "mvtnorm".
##'
##' This test is different from the MAX3 test using in the function
##' \link{max3}. On one hand, the NMAX3 applies to the quantitative
##' traits association studies. However, the MAX3 is used in the
##' case-control association studies. On the other hand, the NMAX3 is
##' based on the nonparametric trend test. However, the MAX3 is based
##' on the Cochran-Armitage trend test.
##'
##' @title The NMAX3 based on the nonparametric trend test in a
##' quantitative trait association study
##' @param y a numeric vector of the observed quantitative trait
##' values in which the \emph{i}th element is the trait value of the
##' \emph{i}th subject.
##' @param g a numeric vector of the observed genotype values (0, 1,
##' or 2 denotes the number of risk alleles) in which the \emph{i}th
##' element is the genotype value of the \emph{i}th subject for a
##' biallelic SNP. \code{g} has the same length as \code{y}.
##' @return A list of \code{test.stat} and
##' \code{p.val}. \code{test.stat} is the observed value of the test
##' statistic and \code{p.val} is the p-value of the test.
##' @author Lin Wang, Wei Zhang, and Qizhai Li.
##' @references W Zhang and Q Li. Nonparametric risk and nonparametric
##' odds in quantitative genetic association studies. \emph{Science
##' Reports (2nd revision)}. 2015.
##' @references B Freidlin, G Zheng, Z Li, and JL Gastwirth. Trend
##' tests for case-control studies of genetic markers: power, sample
##' size and robustness. \emph{Human Heredity}. 2002; 53:146-152.
##' @references WG Cochran. Some methods for strengthening the common
##' chi-square tests. \emph{Biometrics}. 1954; 10:417-451.
##' @references P Armitage. Tests for linear trends in proportions and
##' frequencies. \emph{Biometrics}. 1955; 11:375-386.
##' @examples
##' g <- rbinom(1500, 2, 0.3)
##' y <- 0.5 + 0.25 * g + rgev(1500, 0, 0, 5)
##' nmax3(y, g)
##' @export
nmax3 <- function(y, g)
{
    Z.R <- npt(y, g, 0)[[1]]
    Z.A <- npt(y, g, 0.5)[[1]]
    Z.D <- npt(y, g, 1)[[1]]

    T.max <- max(abs(Z.R), abs(Z.A), abs(Z.D))
    corr.mat <- CorrMatNRTest(y, g)

### calculate the p-value
    pval <- 1-mvtnorm::pmvnorm(lower=-c(T.max,T.max,T.max),upper=c(T.max,T.max,T.max),sigma=corr.mat)[1]

    list(test.stat=T.max, p.val=pval)

}
