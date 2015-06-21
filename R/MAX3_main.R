##' Conduct MAX3 (the maximal value of the three Cochran-Armitage
##' trend tests derived for the recessive, additive, and dominant
##' models) based on the trend tests without the adjustment of the
##' covariates or based on the Wald tests with the adjustment of the
##' covariates to test for the association between a single-nucleotide
##' polymorphism and the binary phenotype.
##'
##' In an association study, the genetic inheritance models
##' (recessive, additive, or dominant) are unknown beforehand. This
##' function can account for the uncertainty of the underlying genetic
##' models and test for the association between a single-nucleotide
##' polymorphism and a binary phenotype with or without correcting for
##' the covariates.
##' @title Maximum Test: maximal value of the three Cochran-Armitage
##' trend tests under the recessive, additive, and dominant models
##' @param y a numeric vector of the observed trait values in which
##' the \emph{i}th element is for the \emph{i}th subject. The elements
##' should be 0 or 1.
##' @param g a numeric vector of the observed genotype values (0, 1,
##' or 2 denotes the number of risk alleles) in which the \emph{i}th
##' element is for the \emph{i}th subject.  The missing value is
##' represented by NA. \code{g} has the same length as \code{y}.
##' @param covariates a numeric matrix for the covariates used in the
##' model. Each column is for one covariate. The default is NULL, that
##' is, there are no covariates to be adjusted for.
##' @param Score.test logical. If TRUE, the score tests are used. One
##' of the \code{Score.test} and the \code{Wald.test} should be FALSE,
##' and the other should be TRUE. The default is TRUE.
##' @param Wald.test logical. If TRUE, the Wald tests are used. One of
##' the \code{Score.test} and the \code{Wald.test} should be FALSE,
##' and the other should be TRUE. The default is FALSE.
##' @param rhombus.formula logical. If TRUE, The p-value for the MAX3
##' is approximated by the rhombus formula. IF FALSE, the 2-fold
##' integration is used to calculate the p-value. The default is
##' FALSE.
##' @return A list of \code{test.stat} and
##' \code{p.val}. \code{test.stat} is the observed value of the test
##' statistic and \code{p.val} is the p-value of the test.
##' @author Lin Wang, Wei Zhang, and Qizhai Li.
##' @references Q Li, G Zheng, Z Li, and K Yu. Efficient approximation
##' of p-value of the maximum of correlated tests, with applications
##' to genome-wide association studies. \emph{Annals of Human
##' Genetics}. 2008; 72(3): 397-406.
##' @examples
##' y <- rep(c(0, 1), 5)
##' g <- sample(c(0, 1, 2), 10, replace = TRUE)
##' max3(y, g, covariates = NULL, Score.test = TRUE, Wald.test = FALSE,
##'        rhombus.formula = FALSE)
##' max3(y, g, covariates = matrix(sample(c(0,1), 20, replace = TRUE), ncol=2),
##'        Score.test = TRUE, Wald.test = FALSE, rhombus.formula = FALSE)
##' @export
max3 <- function(y, g, covariates = NULL, Score.test = TRUE, Wald.test = FALSE, rhombus.formula = FALSE)
{
    dex <- which(!is.na(g))
    g <- g[dex]
    y <- y[dex]

    if (!is.null(covariates))
    {
        covariates <- as.matrix(covariates[dex,])
        N <- length(g)
        g <- matrix(data=g, ncol=1)
        outcome <- rep(y, each=3)
        x.mat <- ChangeX(N, g, covariates=covariates, num.test=3)

        L <- ncol(covariates)
        C.mat <- matrix(data=0, nrow=3, ncol=ncol(x.mat))
        C.mat[1,L+2] <- 1
        C.mat[2,2*L+4] <- 1
        C.mat[3,3*L+6] <- 1

        if (Wald.test)
        {
            temp <- WaldTest(x.mat, outcome, num.test=3, C.mat)
            T.max <- max(abs(temp[[1]]))
            cor.rec.add <- temp[[2]][1,2]
            cor.rec.dom <- temp[[2]][1,3]
            cor.add.dom <- temp[[2]][2,3]
            if (rhombus.formula)
            {
                cor1 <- c(1, cor.rec.add, cor.rec.dom)
                cor2 <- c(cor.rec.add, 1, cor.add.dom)
                cor3 <- c(cor.rec.dom, cor.add.dom, 1)
                cor.matrix <- rbind(cor1, cor2, cor3)
                p.value <- RhombusFormula(cor.matrix, T.max)
            }else
            {
                p.value <- max3Sign(T.max, cor.rec.add, cor.rec.dom, cor.add.dom)
            }
        }else
        {
            id.null <- c(1:(L+1), (L+3):(2*L+3), (2*L+5):(3*L+5))
            temp <- ScoreTest(x.mat, outcome, num.test=3, C.mat, id.null)
            T.max <- max(abs(temp[[1]]))
            cor.rec.add <- temp[[2]][1,2]
            cor.rec.dom <- temp[[2]][1,3]
            cor.add.dom <- temp[[2]][2,3]
            if (rhombus.formula)
            {
                cor1 <- c(1, cor.rec.add, cor.rec.dom)
                cor2 <- c(cor.rec.add, 1, cor.add.dom)
                cor3 <- c(cor.rec.dom, cor.add.dom, 1)
                cor.matrix <- rbind(cor1, cor2, cor3)
                p.value <- RhombusFormula(cor.matrix, T.max)
            }else
            {
                p.value <- max3Sign(T.max, cor.rec.add, cor.rec.dom, cor.add.dom)
            }
        }
    }else
    {
        ri <- c(sum(g==0 & y==1), sum(g==1 & y==1), sum(g==2 & y==1))
        si <- c(sum(g==0 & y==0), sum(g==1 & y==0), sum(g==2 & y==0))

#        T.rec <- TrendTest(0.0, ri=ri, si=si)
#        T.add <- TrendTest(0.5, ri=ri, si=si)
#        T.dom <- TrendTest(1.0, ri=ri, si=si)
        T.rec <- TrendTest(ri=ri, si=si, 0.0)$test.stat
        T.add <- TrendTest(ri=ri, si=si, 0.5)$test.stat
        T.dom <- TrendTest(ri=ri, si=si, 1.0)$test.stat
        T.max <- max(T.rec, T.add, T.dom)

        p <- (ri+si)/(sum(ri+si)+1e-10)
        p0 <- ifelse(p[1]<1e-10, 1e-10, p[1])
        p0 <- ifelse(p[1]>1-1e-10, 1-1e-10, p[1])
        p1 <- ifelse(p[2]<1e-10, 1e-10, p[2])
        p1 <- ifelse(p[2]>1-1e-10, 1-1e-10, p[2])
        p2 <- ifelse(p[3]<1e-10, 1e-10, p[3])
        p2 <- ifelse(p[3]>1-1e-10, 1-1e-10, p[3])

        cor.rec.add <- p2*(2*p0+p1)/sqrt(p2*(1-p2))/sqrt(p0*(p1+2*p2)+p2*(p1+2*p0))
        cor.rec.dom <- p0*p2/sqrt(p0*(1-p0))/sqrt(p2*(1-p2))
        cor.add.dom <- p0*(p1+2*p2)/sqrt(p0*(1-p0))/sqrt(p0*(p1+2*p2)+p2*(p1+2*p0))
        cor.matrix <- matrix(data=c(1, cor.rec.add, cor.rec.dom, cor.rec.add, 1, cor.add.dom, cor.rec.dom, cor.add.dom, 1), nrow=3, ncol=3, byrow=TRUE)
        #1-pmvnorm(lower=-c(T.max,T.max,T.max),upper=c(T.max,T.max,T.max),sigma=cor.matrix)

        if (rhombus.formula)
        {
            p.value <- RhombusFormula(cor.matrix, T.max)
        }else
        {
            p.value <- max3Sign(T.max, cor.rec.add, cor.rec.dom, cor.add.dom)
        }
    }

    pv <- min(1,p.value)
    list(test.stat=T.max, p.val=pv)
}
