% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/es_from_ANOVA_means.R
\name{es_from_means_sd_pooled}
\alias{es_from_means_sd_pooled}
\title{Convert means of two groups and the pooled standard deviation into several effect size measures}
\usage{
es_from_means_sd_pooled(
  mean_exp,
  mean_nexp,
  mean_sd_pooled,
  n_exp,
  n_nexp,
  smd_to_cor = "viechtbauer",
  reverse_means
)
}
\arguments{
\item{mean_exp}{mean of participants in the experimental/exposed group.}

\item{mean_nexp}{mean of participants in the non-experimental/non-exposed group.}

\item{mean_sd_pooled}{pooled standard deviation across both groups.}

\item{n_exp}{number of participants in the experimental/exposed group.}

\item{n_nexp}{number of participants in the non-experimental/non-exposed group.}

\item{smd_to_cor}{formula used to convert the \code{cohen_d} value into a coefficient correlation (see details).}

\item{reverse_means}{a logical value indicating whether the direction of the generated effect sizes should be flipped.}
}
\value{
This function estimates and converts between several effect size measures.

\tabular{ll}{
\code{natural effect size measure} \tab MD + D + G\cr
\tab \cr
\code{converted effect size measure} \tab OR + R + Z \cr
\tab \cr
\code{required input data} \tab See 'Section 9. Means and dispersion (crude)'\cr
\tab https://metaconvert.org/html/input.html\cr
\tab \cr
}
}
\description{
Convert means of two groups and the pooled standard deviation into several effect size measures
}
\details{
This function first computes a Cohen's d (D), Hedges' g (G) and mean difference (MD)
from the means of two independent groups and the pooled standard deviation across the groups.
Odds ratio (OR) and correlation coefficients (R/Z) are then converted from the Cohen's d.

\strong{To estimate a mean difference}  (formulas 12.1-12.6 in Cooper):
\deqn{md = mean\_exp - mean\_nexp}
\deqn{md\_se = \sqrt{\frac{n_exp+n_nexp}{n_exp*n_nexp} * mean_sd_pooled^2}}
\deqn{md\_ci\_lo = md - md\_se * qt(.975, df = n\_exp + n\_nexp - 2)}
\deqn{md\_ci\_up = md + md\_se * qt(.975, df = n\_exp + n\_nexp - 2)}

\strong{To estimate a Cohen's d} the following formulas are used (formulas 12.10-12.18 in Cooper):
\deqn{cohen\_d =  \frac{mean\_exp - mean\_nexp}{means\_sd\_pooled}}
\deqn{cohen\_d\_se = \frac{(n\_exp+n\_nexp)}{n\_exp*n\_nexp} + \frac{cohen\_d^2}{2(n\_exp+n\_nexp)}}
\deqn{cohen\_d\_ci\_lo = cohen\_d - cohen\_d\_se * qt(.975, df = n\_exp + n\_nexp - 2)}
\deqn{cohen\_d\_ci\_up = cohen\_d + cohen\_d\_se * qt(.975, df = n\_exp + n\_nexp - 2)}

\strong{To estimate other effect size measures},
calculations of the \code{\link{es_from_cohen_d}()} are applied.
}
\examples{
es_from_means_sd_pooled(
  n_exp = 55, n_nexp = 55,
  mean_exp = 2.3, mean_nexp = 1.9,
  mean_sd_pooled = 0.9
)
}
\references{
Cooper, H., Hedges, L.V., & Valentine, J.C. (Eds.). (2019). The handbook of research synthesis and meta-analysis. Russell Sage Foundation.
}
