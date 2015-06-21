##' Conduct the single-marker test in an association study to test for
##' the association between the genotype at a biallelic marker and a
##' trait.
##'
##' Single-marker analysis is a core in many gene-based or
##' pathway-based procedures, such as the truncated p-value
##' combination and the minimal p-value.
##' @title Single-marker test
##' @param y a numeric vector of the observed trait values in which
##' the \emph{i}th element is for the \emph{i}th subject. The elements
##' could be discrete (0 or 1) or continuous. The missing value is
##' represented by NA.
##' @param g a numeric vector of the observed genotype values (0, 1,
##' or 2 denotes the number of risk alleles) in which the \emph{i}th
##' element is for the \emph{i}th subject.  The missing value is
##' represented by NA. \code{g} has the same length as \code{y}.
##' @param covariates an optional data frame, list or environment
##' containing the covariates used in the model. The default is NULL,
##' that is, there are no covariates.
##' @param min.count a critical value to decide which method is used
##' to calculate the p-value when the trait is discrete and covariates
##' = NULL. If the minimum number of the elements given a specific
##' trait value and a specific genotype value is less than
##' \code{min.count}, the Fisher exact test is adopted; otherwise, the
##' Wald test is adopted.
##' @param missing.rate the highest missing rate of the genotype
##' values that this function can tolerate.
##' @param y.continuous logical. If TRUE, \code{y} is continuous;
##' otherwise, \code{y} is discrete.
##' @return \code{smt} returns a list.
##'
##' If y is continuous, the list contains the following components:
##' \tabular{llll}{
##' n0 \tab \tab \tab the number of the subjects used for the calculation with the genotype 0.\cr
##' n1 \tab \tab \tab the number of the subjects used for the calculation with the genotype 1.\cr
##' n2 \tab \tab \tab the number of the subjects used for the calculation with the genotype 2.\cr
##' p.val \tab \tab \tab the p-value of the single-marker test.
##' }
##' If y is discrete, the list contains the following components:
##' \tabular{llll}{
##' r0 \tab \tab \tab the number of the subjects used for the calculation with the trait value 1 and the genotype 0.\cr
##' r1 \tab \tab \tab the number of the subjects used for the calculation with the trait value 1 and the genotype 1.\cr
##' r2 \tab \tab \tab the number of the subjects used for the calculation with the trait value 1 and the genotype 2.\cr
##' r.miss \tab \tab \tab the number of the missing genotype values with the trait value 1.\cr
##' s0 \tab \tab \tab the number of the subjects used for the calculation with the trait value 0 and the genotype 0.\cr
##' s1 \tab \tab \tab the number of the subjects used for the calculation with the trait value 0 and the genotype 1.\cr
##' s2 \tab \tab \tab the number of the subjects used for the calculation with the trait value 0 and the genotype 2.\cr
##' s.miss \tab \tab \tab the number of the missing genotype values with the trait value 0.\cr
##' n.miss \tab \tab \tab the total number of the missing genotype values. \cr
##' p.val \tab \tab \tab the p-value of the single-marker test.
##' }
##' @author Lin Wang, Wei Zhang, and Qizhai Li.
##' @examples
##' y <- rep(c(0, 1), 25)
##' g <- sample(c(0, 1, 2), 50, replace = TRUE)
##' smt(y, g, covariates = NULL, min.count=5,
##'         missing.rate=0.20, y.continuous = FALSE)
##' @export
smt <- function(y, g, covariates=NULL, min.count=5, missing.rate=0.20, y.continuous=FALSE)
{
    # g is the genotype vector taking 0, 1 and 2
    # y is the outcome

    dex <- !is.na(g) & !is.na(y)
    g.valid <- g[dex]
    y.valid <- y[dex]

    # y is continuous
    if (y.continuous)
    {
        n0 <- sum(g.valid==0)
        n1 <- sum(g.valid==1)
        n2 <- sum(g.valid==2)

        if(is.null(covariates))
        {
            temp <- lm(y~g)
        }
        else
        {
            temp <- lm(y~.+g, data=covariates)
        }
        a.1 <- summary(temp)$coefficients
        t.1 <- dim(a.1)
        p.value <- a.1[t.1[1], t.1[2]]

        res <- list(n0=n0, n1=n1, n2=n2, p.val=p.value)
        return(res)
    }
    else
    {
        r0 <- sum(g.valid==0 & y.valid==1)
        r1 <- sum(g.valid==1 & y.valid==1)
        r2 <- sum(g.valid==2 & y.valid==1)

        s0 <- sum(g.valid==0 & y.valid==0)
        s1 <- sum(g.valid==1 & y.valid==0)
        s2 <- sum(g.valid==2 & y.valid==0)

        addr <- !is.na(y)
        u <- g[addr]
        v <- y[addr]
        r.miss <- sum(is.na(u[v==1]))
        s.miss <- sum(is.na(u[v==0]))

        n.miss <- sum(is.na(g))

        # missing rate
        m.r <- sum(is.na(g))/length(g)
        if (m.r>=missing.rate)
        {
            res <- c(r0=r0, r1=r1, r2=r2, r.miss=r.miss, s0=s0, s1=s1, s2=s2, s.miss=s.miss, n.miss=n.miss, p.val=-9999)
            return(res)
        }
        else
        {
            if (min(CalExpect(matrix(c(r0,r1,r2,s0,s1,s2),nrow=2,byrow=TRUE)))>=min.count)
            {
                if (is.null(covariates))
                {
                    temp <- glm(y~g, family=binomial(link="logit"))
                }
                else
                {
                    temp <- glm(y~.+g, family=binomial(link="logit"), data=covariates)
                }
                a.1 <- summary(temp)$coefficients
                t.1 <- dim(a.1)
                p.value <- a.1[t.1[1], t.1[2]]
            }
            else
            {
                if (s0 >= s2)
                {
                    H <- matrix(c(r0,r1+r2,s0,s1+s2),ncol=2, byrow=TRUE)
                    if (min(CalExpect(H))<min.count)
                    {
                        p.value <- fisher.test(H)$p.value
                    }
                    else
                    {
                        g[g==2] <- 1
                        if (is.null(covariates))
                        {
                            temp <- glm(y~g, family=binomial(link="logit"))
                        }
                        else
                        {
                            temp <- glm(y~.+g, family=binomial(link="logit"), data=covariates)
                        }
                        a.1 <- summary(temp)$coefficients
                        t.1 <- dim(a.1)
                        p.value <- a.1[t.1[1], t.1[2]]
                    }
                }
                else
                {
                    H <- matrix(c(r0+r1,r2,s0+s1,s2),ncol=2,byrow=TRUE)
                    if (min(CalExpect(H))<min.count)
                    {
                        p.value <- fisher.test(H)$p.value
                    }
                    else
                    {
                        g[g==1] <- 0
                        g[g==2] <- 1
                        if (is.null(covariates))
                        {
                            temp <- glm(y~g, family=binomial(link="logit"))
                        }
                        else
                        {
                            temp <- glm(y~.+g, family=binomial(link="logit"), data=covariates)
                        }
                        a.1 <- summary(temp)$coefficients
                        t.1 <- dim(a.1)
                        p.value <- a.1[t.1[1], t.1[2]]
                    }
                }
            }
            res <- list(r0=r0, r1=r1, r2=r2, r.miss=r.miss, s0=s0, s1=s1, s2=s2, s.miss=s.miss, n.miss=n.miss, p.val=p.value)
            return(res)
        }
    }
}
