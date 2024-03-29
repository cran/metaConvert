% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/es_from_ANOVA_med_range_quart.R
\name{es_from_med_quarts}
\alias{es_from_med_quarts}
\title{Convert median and interquartile range of two independent groups into several effect size measures}
\usage{
es_from_med_quarts(
  q1_exp,
  med_exp,
  q3_exp,
  n_exp,
  q1_nexp,
  med_nexp,
  q3_nexp,
  n_nexp,
  smd_to_cor = "viechtbauer",
  reverse_med
)
}
\arguments{
\item{q1_exp}{first quartile of the experimental/exposed group.}

\item{med_exp}{median value of the experimental/exposed group.}

\item{q3_exp}{third quartile of the experimental/exposed group.}

\item{n_exp}{number of participants in the experimental/exposed group.}

\item{q1_nexp}{first quartile of the non-experimental/non-exposed group.}

\item{med_nexp}{median value of the non-experimental/non-exposed group.}

\item{q3_nexp}{third quartile of the non-experimental/non-exposed group.}

\item{n_nexp}{number of participants in the non-experimental/non-exposed group.}

\item{smd_to_cor}{formula used to convert the generated \code{cohen_d} value into a coefficient correlation (see details).}

\item{reverse_med}{a logical value indicating whether the direction of generated effect sizes should be flipped.}
}
\value{
This function estimates and converts between several effect size measures.
}
\description{
Convert median and interquartile range of two independent groups into several effect size measures
}
\details{
This function first converts a Cohen's d (D), Hedges' g (G) and mean difference (MD)
from the medians and interquartile ranges of two independent groups.
Odds ratio (OR) and correlation coefficients (R/Z) are then converted from the Cohen's d.

\strong{This function recreates means+SD of the two groups}  (Wan et al., 2014):
\deqn{mean\_exp = \frac{q1\_exp + med\_exp + q3\_exp}{3}}
\deqn{mean\_nexp = \frac{q1\_nexp + med\_nexp + q3\_nexp}{3}}
\deqn{mean\_sd\_exp = \frac{q3\_exp - q1\_exp}{2*qnorm(\frac{0.75*n\_exp - 0.125}{n\_exp+0.25})}}
\deqn{mean\_sd\_nexp = \frac{q3\_nexp - q1\_nexp}{2*qnorm(\frac{0.75*n\_nexp - 0.125}{n\_nexp+0.25})}}

Note that if the group sample size is inferior to 50, a correction is applied to estimate the standard deviation.

\strong{From these means+SD, the function computes MD, D and G} using formulas
described in \code{\link{es_from_means_sd}()}.

\strong{To estimate other effect size measures},
calculations of the \code{\link{es_from_cohen_d}()} are applied.
}
\examples{
es_from_med_quarts(
  q1_exp = 1335, med_exp = 1400,
  q3_exp = 1765, n_exp = 40,
  q1_nexp = 1481, med_nexp = 1625,
  q3_nexp = 1800, n_nexp = 40
)
}
\references{
Wan, X., Wang, W., Liu, J. et al. Estimating the sample mean and standard deviation from the sample size, median, range and/or interquartile range. BMC Med Res Methodol 14, 135 (2014). https://doi.org/10.1186/1471-2288-14-135
}
