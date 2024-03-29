% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/es_from_PHI_CHISQ.R
\name{es_from_chisq}
\alias{es_from_chisq}
\title{Convert a chi-square value to several effect size measures}
\usage{
es_from_chisq(
  chisq,
  n_sample,
  n_cases,
  n_exp,
  yates_chisq = FALSE,
  reverse_chisq
)
}
\arguments{
\item{chisq}{value of the chi-squared}

\item{n_sample}{total number of participants in the sample}

\item{n_cases}{total number of cases/events}

\item{n_exp}{total number of participants in the exposed group}

\item{yates_chisq}{a logical value indicating whether the Chi square has been performed using Yate's correction for continuity.}

\item{reverse_chisq}{a logical value indicating whether the direction of generated effect sizes should be flipped.}
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
Convert a chi-square value to several effect size measures
}
\details{
This function converts a chi-square value (with one degree of freedom)
into a phi coefficient (Lipsey et al. 2001):
\deqn{phi = chisq * \sqrt{\frac{chisq}{n\_sample}}}.

Note that if \code{yates_chisq = "TRUE"}, a small correction is added.

Then, the phi coefficient is converted to other effect size measures (see \code{\link{es_from_phi}}).
}
\examples{
es_from_chisq(chisq = 4.21, n_sample = 78, n_cases = 51, n_exp = 50)
}
\references{
Lipsey, M. W., & Wilson, D. B. (2001). Practical meta-analysis. Sage Publications, Inc.
}
