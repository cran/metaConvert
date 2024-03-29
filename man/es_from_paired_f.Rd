% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/es_from_PAIRED_STATISTICS.R
\name{es_from_paired_f}
\alias{es_from_paired_f}
\title{Convert two paired ANOVA f value of two independent groups into several effect size measures}
\usage{
es_from_paired_f(
  paired_f_exp,
  paired_f_nexp,
  n_exp,
  n_nexp,
  r_pre_post_exp,
  r_pre_post_nexp,
  smd_to_cor = "viechtbauer",
  reverse_paired_f
)
}
\arguments{
\item{paired_f_exp}{Paired ANOVA F value of the experimental/exposed group.}

\item{paired_f_nexp}{Paired ANOVA F value of the non-experimental/non-exposed group.}

\item{n_exp}{number of participants in the experimental/exposed group.}

\item{n_nexp}{number of participants in the non-experimental/non-exposed group.}

\item{r_pre_post_exp}{pre-post correlation in the experimental/exposed group}

\item{r_pre_post_nexp}{pre-post correlation in the non-experimental/non-exposed group}

\item{smd_to_cor}{formula used to convert the \code{cohen_d} value into a coefficient correlation (see details).}

\item{reverse_paired_f}{a logical value indicating whether the direction of generated effect sizes should be flipped.}
}
\value{
This function estimates and converts between several effect size measures.

\tabular{ll}{
\code{natural effect size measure} \tab D + G\cr
\tab \cr
\code{converted effect size measure} \tab OR + R + Z \cr
\tab \cr
\code{required input data} \tab See 'Section 16. Paired: Paired F- or t-test'\cr
\tab https://metaconvert.org/html/input.html\cr
\tab \cr
}
}
\description{
Convert two paired ANOVA f value of two independent groups into several effect size measures
}
\details{
This function converts the paired F-test obtained from two independent groups value into a Cohen's d (D) and Hedges' g (G) (table 12.2 in Cooper).
Odds ratio (OR) and correlation coefficients (R/Z) are then converted from the Cohen's d.

\strong{To estimate the Cohen's d,} the following formulas are used (Cooper et al., 2019):
This function converts a Student's t-test value into a Cohen's d (table 12.2 in Cooper).
\deqn{paired\_t\_exp = \sqrt{paired\_f\_exp}}
\deqn{paired\_t\_nexp = \sqrt{paired\_f\_nexp}}

\strong{To estimate other effect size measures},
calculations of the \code{\link{es_from_paired_t}()} are applied.
}
\examples{
es_from_paired_f(paired_f_exp = 2.1, paired_f_nexp = 4.2, n_exp = 20, n_nexp = 22)
}
\references{
Cooper, H., Hedges, L.V., & Valentine, J.C. (Eds.). (2019). The handbook of research synthesis and meta-analysis. Russell Sage Foundation.
}
