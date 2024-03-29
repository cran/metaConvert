% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/es_from_PHI_CHISQ.R
\name{es_from_chisq_pval}
\alias{es_from_chisq_pval}
\title{Convert a p-value of a chi-square to several effect size measures}
\usage{
es_from_chisq_pval(
  chisq_pval,
  n_sample,
  n_cases,
  n_exp,
  yates_chisq = FALSE,
  reverse_chisq_pval
)
}
\arguments{
\item{chisq_pval}{p-value of a chi-square coefficient}

\item{n_sample}{total number of participants in the sample}

\item{n_cases}{total number of cases/events}

\item{n_exp}{total number of participants in the exposed group}

\item{yates_chisq}{a logical value indicating whether the Chi square has been performed using Yate's correction for continuity.}

\item{reverse_chisq_pval}{a logical value indicating whether the direction of generated effect sizes should be flipped.}
}
\value{
This function estimates and converts between several effect size measures.

\tabular{ll}{
\code{natural effect size measure} \tab OR + RR + NNT\cr
\tab \cr
\code{converted effect size measure} \tab D + G + R + Z \cr
\tab \cr
\code{required input data} \tab See 'Section 8. Phi or chi-square'\cr
\tab https://metaconvert.org/html/input.html\cr
\tab \cr
}
}
\description{
Convert a p-value of a chi-square to several effect size measures
}
\details{
This function converts a chi-square value (with one degree of freedom)
into a chi-square coefficient (Section 3.12 in Lipsey et al., 2001):
\deqn{chisq = qchisq(chisq\_pval, df = 1, lower.tail = FALSE)}

Note that if \code{yates_chisq = "TRUE"}, a small correction is added.

Then, the chisq coefficient is converted to other effect size measures (see \code{\link{es_from_chisq}}).
}
\examples{
es_from_chisq_pval(chisq_pval = 0.2, n_sample = 42, n_exp = 25, n_cases = 13)
}
\references{
Lipsey, M. W., & Wilson, D. B. (2001). Practical meta-analysis. Sage Publications, Inc.
}
