% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/es_from_stand_SMD.R
\name{es_from_hedges_g}
\alias{es_from_hedges_g}
\title{Convert a Hedges' g value to other effect size measures (G, OR, COR)}
\usage{
es_from_hedges_g(
  hedges_g,
  n_exp,
  n_nexp,
  smd_to_cor = "viechtbauer",
  reverse_g
)
}
\arguments{
\item{hedges_g}{Hedges' g value}

\item{n_exp}{number of participants in the experimental/exposed group.}

\item{n_nexp}{number of participants in the non-experimental/non-exposed group.}

\item{smd_to_cor}{formula used to convert the \code{hedges_g} value into a coefficient correlation (see details).}

\item{reverse_g}{a logical value indicating whether the direction of the \code{hedges_g} value should be flipped.}
}
\value{
This function estimates and converts between several effect size measures.

\tabular{ll}{
\code{natural effect size measure} \tab D + G\cr
\tab \cr
\code{converted effect size measure} \tab OR + R + Z\cr
\tab \cr
\code{required input data} \tab See 'Section 1. Cohen's d or Hedges' g'\cr
\tab https://metaconvert.org/html/input.html\cr
\tab \cr
}
}
\description{
Convert a Hedges' g value to other effect size measures (G, OR, COR)
}
\details{
This function estimates the standard error of the Hedges' g and the Cohen's d (D).
Odds ratio (OR) and correlation coefficients (R/Z) are then converted from the Cohen's d.

\strong{To estimate standard error of Hedges'g}, the following formula is used (Hedges, 1981):
\deqn{df = n\_exp + n\_nexp - 2}
\deqn{hedges\_g\_se = \sqrt{cohen\_d\_se^2 * J^2}}
\deqn{hedges\_g\_ci\_lo = hedges\_g - hedges\_g\_se * qt(.975, df = n\_exp+n\_nexp-2)}
\deqn{hedges\_g\_ci\_up = hedges\_g + hedges\_g\_se * qt(.975, df = n\_exp+n\_nexp-2)}

\strong{To estimate the Cohen's d value}, the following formula is used (Hedges, 1981):
\deqn{J = exp(\log_{gamma}(\frac{df}{2}) - 0.5 * \log(\frac{df}{2}) - \log_{gamma}(\frac{df - 1}{2}))}
\deqn{cohen\_d = \frac{hedges\_g}{J}}
\deqn{cohen\_d\_se = \sqrt{(\frac{n\_exp+n\_nexp}{n\_exp*n\_nexp} + \frac{cohen\_d^2}{2*(n\_exp+n\_nexp)})}}

\strong{To estimate other effect size measures},
calculations of the \code{\link{es_from_cohen_d}()} are applied.
}
\examples{
es_from_hedges_g(hedges_g = 0.243, n_exp = 20, n_nexp = 20)
}
\references{
Hedges LV (1981): Distribution theory for Glass’s estimator of effect size and related estimators. Journal of Educational and Behavioral Statistics, 6, 107–28
}
