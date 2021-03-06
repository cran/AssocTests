% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SMT_main.R
\name{smt}
\alias{smt}
\title{Single-marker test}
\usage{
smt(
  y,
  g,
  covariates = NULL,
  min.count = 5,
  missing.rate = 0.2,
  y.continuous = FALSE
)
}
\arguments{
\item{y}{a numeric vector of the observed trait values in which
the \emph{i}th element is for the \emph{i}th subject. The elements
could be discrete (\code{0} or \code{1}) or continuous. The missing value is
represented by \code{NA}.}

\item{g}{a numeric vector of the observed genotype values (\code{0}, \code{1},
or \code{2} denotes the number of risk alleles) in which the \emph{i}th
element is for the \emph{i}th subject. The missing value is
represented by \code{NA}. \code{g} has the same length as \code{y}.}

\item{covariates}{an optional data frame, list or environment
containing the covariates used in the model. The default is \code{NULL},
that is, there are no covariates.}

\item{min.count}{a critical value to decide which method is used
to calculate the p-value when the trait is discrete and \code{covariates
= NULL}. If the minimum number of the elements given a specific
trait value and a specific genotype value is less than
\code{min.count}, the Fisher's exact test is adopted; otherwise, the
Wald test is adopted. The default is \code{5}.}

\item{missing.rate}{the highest missing value rate of the genotype
values that this function can tolerate. The default is \code{0.2}.}

\item{y.continuous}{logical. If \code{TRUE}, \code{y} is continuous;
otherwise, \code{y} is discrete. The default is \code{FALSE}.}
}
\value{
\code{smt} returns a list with class "\code{htest}".

If y is continuous, the list contains the following components:
\tabular{llll}{
\code{statistic} \tab \tab \tab \cr
\tab \tab \tab the observed value of the test statistic.\cr
\code{p.value} \tab \tab \tab \cr
\tab \tab \tab the p-value for the test.\cr
\code{alternative} \tab \tab \tab \cr
\tab \tab \tab a character string describing the alternative hypothesis.\cr
\code{method} \tab \tab \tab \cr
\tab \tab \tab a character string indicating the type of test performed.\cr
\code{data.name} \tab \tab \tab \cr
\tab \tab \tab a character string giving the names of the data. \cr
\code{sample.size} \tab \tab \tab \cr
\tab \tab \tab a vector giving the numbers of the subjects with the genotypes \code{0}, \code{1}, and \code{2} (\code{n0}, \cr
\tab \tab \tab \code{n1}, and \code{n2}, respectively).
}
If y is discrete, the list contains the following components:
\tabular{llll}{
\code{statistic} \tab \tab \tab \cr
\tab \tab \tab the observed value of the test statistic.\cr
\code{p.value} \tab \tab \tab \cr
\tab \tab \tab the p-value for the test.\cr
\code{alternative} \tab \tab \tab \cr
\tab \tab \tab a character string describing the alternative hypothesis.\cr
\code{method} \tab \tab \tab \cr
\tab \tab \tab a character string indicating the type of test performed.\cr
\code{data.name} \tab \tab \tab \cr
\tab \tab \tab a character string giving the names of the data. \cr
\code{sample.size} \tab \tab \tab \cr
\tab \tab \tab a vector giving \cr
\tab \tab \tab the number of subjects with the trait value \code{1} and the genotype \code{0} (\code{r0}), \cr
\tab \tab \tab the number of subjects with the trait value \code{1} and the genotype \code{1} (\code{r1}), \cr
\tab \tab \tab the number of subjects with the trait value \code{1} and the genotype \code{2} (\code{r2}), \cr
\tab \tab \tab the number of subjects with the trait value \code{0} and the genotype \code{0} (\code{s0}), \cr
\tab \tab \tab the number of subjects with the trait value \code{0} and the genotype \code{1} (\code{s1}), \cr
\tab \tab \tab and the number of subjects with the trait value \code{0} and the genotype \code{2} (\code{s2}).\cr
\code{bad.obs} \tab \tab \tab \cr
\tab \tab \tab a vector giving the number of missing genotype values with the trait value \code{1} \cr
\tab \tab \tab (\code{r.miss}), the number of missing genotype values with the trait value \code{0} \cr
\tab \tab \tab (\code{s.miss}), and the total number of the missing genotype values (\code{n.miss}).
}
}
\description{
Conduct the single-marker test in an association study to test for
the association between the genotype at a biallelic marker and a
trait.
}
\details{
Single-marker analysis is a core in many gene-based or
pathway-based procedures, such as the truncated p-value
combination and the minimal p-value.
}
\examples{
y <- rep(c(0, 1), 25)
g <- sample(c(0, 1, 2), 50, replace = TRUE)
smt(y, g, covariates = NULL, min.count=5,
        missing.rate=0.20, y.continuous = FALSE)
}
\references{
Lin Wang, Wei Zhang, and Qizhai Li. AssocTests: An R Package 
for Genetic Association Studies. \emph{Journal of Statistical Software}. 
2020; 94(5): 1-26.
}
\author{
Lin Wang, Wei Zhang, and Qizhai Li.
}
