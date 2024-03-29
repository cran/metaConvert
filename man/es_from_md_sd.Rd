% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/es_from_stand_MD.R
\name{es_from_md_sd}
\alias{es_from_md_sd}
\title{Convert a mean difference between two independent groups and standard deviation into several effect size measures}
\usage{
es_from_md_sd(md, md_sd, n_exp, n_nexp, smd_to_cor = "viechtbauer", reverse_md)
}
\arguments{
\item{md}{mean difference between two independent groups}

\item{md_sd}{standard deviation of the mean difference}

\item{n_exp}{number of participants in the experimental/exposed group.}

\item{n_nexp}{number of participants in the non-experimental/non-exposed group.}

\item{smd_to_cor}{formula used to convert the \code{cohen_d} value into a coefficient correlation (see details).}

\item{reverse_md}{a logical value indicating whether the direction of generated effect sizes should be flipped.}
}
\value{
This function estimates and converts between several effect size measures.

\tabular{ll}{
\code{natural effect size measure} \tab MD + D + G\cr
\tab \cr
\code{converted effect size measure} \tab OR + R + Z\cr
\tab \cr
\code{required input data} \tab See 'Section 10. Mean difference and dispersion (crude)'\cr
\tab https://metaconvert.org/html/input.html\cr
\tab \cr
}
}
\description{
Convert a mean difference between two independent groups and standard deviation into several effect size measures
}
\details{
This function converts the mean difference and 95\% CI into a Cohen's d (D) and Hedges' g (G).
Odds ratio (OR) and correlation coefficients (R/Z) are then converted from the Cohen's d.

\strong{The formula used to obtain the Cohen's d} is:
\deqn{d = \frac{md}{md\_sd}}
Note that this formula is perfectly accurate only if the \code{md_sd} has been estimated by assuming that the variance of the two groups is equal.

\strong{To estimate other effect size measures},
calculations of the \code{\link{es_from_cohen_d}()} are applied.
}
\examples{
es_from_md_sd(md = 4, md_sd = 2, n_exp = 20, n_nexp = 22)
}
