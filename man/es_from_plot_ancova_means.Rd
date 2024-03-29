% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/es_from_PLOTS.R
\name{es_from_plot_ancova_means}
\alias{es_from_plot_ancova_means}
\title{Converts the means and bounds of an error bar (generally extracted from a plot) into four effect measures (SMD, MD, OR, COR)}
\usage{
es_from_plot_ancova_means(
  n_exp,
  n_nexp,
  plot_ancova_mean_exp,
  plot_ancova_mean_nexp,
  plot_ancova_mean_sd_lo_exp,
  plot_ancova_mean_sd_lo_nexp,
  plot_ancova_mean_sd_up_exp,
  plot_ancova_mean_sd_up_nexp,
  plot_ancova_mean_se_lo_exp,
  plot_ancova_mean_se_lo_nexp,
  plot_ancova_mean_se_up_exp,
  plot_ancova_mean_se_up_nexp,
  plot_ancova_mean_ci_lo_exp,
  plot_ancova_mean_ci_lo_nexp,
  plot_ancova_mean_ci_up_exp,
  plot_ancova_mean_ci_up_nexp,
  cov_outcome_r,
  n_cov_ancova,
  smd_to_cor = "viechtbauer",
  reverse_plot_ancova_means
)
}
\arguments{
\item{n_exp}{number of participants in the experimental/exposed group.}

\item{n_nexp}{number of participants in the non-experimental/non-exposed group.}

\item{plot_ancova_mean_exp}{ancova_mean of participants in the experimental/exposed group (extracted from a plot).}

\item{plot_ancova_mean_nexp}{ancova_mean of participants in the non-experimental/non-exposed group (extracted from a plot).}

\item{plot_ancova_mean_sd_lo_exp}{lower bound of an error bar depicting -1 SD from the ancova_mean of the experimental/exposed group (extracted from a plot).}

\item{plot_ancova_mean_sd_lo_nexp}{lower bound of an error bar depicting -1 SD from the ancova_mean of the non-experimental/non-exposed group (extracted from a plot).}

\item{plot_ancova_mean_sd_up_exp}{upper bound of an error bar depicting +1 SD from the ancova_mean of the experimental/exposed group (extracted from a plot).}

\item{plot_ancova_mean_sd_up_nexp}{upper bound of an error bar depicting +1 SD from the ancova_mean of the non-experimental/non-exposed group (extracted from a plot).}

\item{plot_ancova_mean_se_lo_exp}{lower bound of an error bar depicting -1 SE from the ancova_mean of the experimental/exposed group (extracted from a plot).}

\item{plot_ancova_mean_se_lo_nexp}{lower bound of an error bar depicting -1 SE from the ancova_mean of the non-experimental/non-exposed group (extracted from a plot).}

\item{plot_ancova_mean_se_up_exp}{upper bound of an error bar depicting +1 SE from the ancova_mean of the experimental/exposed group (extracted from a plot).}

\item{plot_ancova_mean_se_up_nexp}{upper bound of an error bar depicting +1 SE from the ancova_mean of the non-experimental/non-exposed group (extracted from a plot).}

\item{plot_ancova_mean_ci_lo_exp}{lower bound of an error bar depicting the 95\% CI of the ancova_mean of the experimental/exposed group (extracted from a plot).}

\item{plot_ancova_mean_ci_lo_nexp}{lower bound of an error bar depicting the 95\% CI of the ancova_mean of the non-experimental/non-exposed group (extracted from a plot).}

\item{plot_ancova_mean_ci_up_exp}{upper bound of an error bar depicting the 95\% CI of the ancova_mean of the experimental/exposed group (extracted from a plot).}

\item{plot_ancova_mean_ci_up_nexp}{upper bound of an error bar depicting the 95\% CI of the ancova_mean of the non-experimental/non-exposed group (extracted from a plot).}

\item{cov_outcome_r}{correlation between the outcome and covariate (multiple correlation when multiple covariates are included in the ANCOVA model).}

\item{n_cov_ancova}{number of covariates in the ANCOVA model.}

\item{smd_to_cor}{formula used to convert the \code{cohen_d} value into a coefficient correlation (see details).}

\item{reverse_plot_ancova_means}{a logical value indicating whether the direction of the generated effect sizes should be flipped.}
}
\value{
This function estimates and converts between several effect size measures.

\tabular{ll}{
\code{natural effect size measure} \tab MD + D + G\cr
\tab \cr
\code{converted effect size measure} \tab OR + R + Z \cr
\tab \cr
\code{required input data} \tab See 'Section 22. From plot: adjusted means and dispersion (adjusted)'\cr
\tab https://metaconvert.org/html/input.html\cr
\tab \cr
}
}
\description{
Converts the means and bounds of an error bar (generally extracted from a plot) into four effect measures (SMD, MD, OR, COR)
}
\details{
This function uses the bounds of an error bar of a mean obtained from a plot into a standard deviation.
Then, a mean difference (MD), Cohen's d (D), and Hedges' g (G) are estimated.
Odds ratio (OR), risk ratio (RR) and correlation coefficients (R/Z) are converted from the Cohen's d value.

\strong{To convert the bound of an error bar into a standard deviation,}
this function always prioritizes information from the \code{plot_ancova_mean_sd_*} arguments,
then those from the \code{plot_ancova_mean_se_*} arguments,
then those from the \code{plot_ancova_mean_ci_*} arguments.
\enumerate{
\item If the bounds of the standard deviations are provided, the following formulas are used:
\deqn{ancova\_mean\_sd\_lo\_exp = plot\_ancova\_mean\_exp - plot\_ancova\_mean\_sd\_lo\_exp}
\deqn{ancova\_mean\_sd\_up\_exp = plot\_ancova\_mean\_sd\_up\_exp - plot\_ancova\_mean\_exp}
\deqn{ancova\_mean\_sd\_exp = \frac{ancova\_mean\_sd\_lo\_exp + ancova\_mean\_sd\_up\_exp}{2}}
}

\deqn{mean\_sd\_lo\_nexp = plot\_ancova\_mean\_nexp - plot\_ancova\_mean\_sd\_lo\_nexp}
\deqn{mean\_sd\_up\_nexp = plot\_ancova\_mean\_sd\_up\_nexp - plot\_ancova\_mean\_nexp}
\deqn{mean\_sd\_nexp = \frac{mean\_sd\_lo\_nexp + mean\_sd\_up\_nexp}{2}}

Then, calculations of the \code{\link{es_from_ancova_means_sd}} are used.
\enumerate{
\item If the bounds of the standard errors are provided, the following formulas are used:
\deqn{ancova\_mean\_se\_lo\_exp = plot\_ancova\_mean\_exp - plot\_ancova\_mean\_se\_lo\_exp}
\deqn{ancova\_mean\_se\_up\_exp = plot\_ancova\_mean\_se\_up\_exp - plot\_ancova\_mean\_exp}
\deqn{ancova\_mean\_se\_exp = \frac{ancova\_mean\_se\_lo\_exp + ancova\_mean\_se\_up\_exp}{2}}
}

\deqn{mean\_se\_lo\_nexp = plot\_ancova\_mean\_nexp - plot\_ancova\_mean\_se\_lo\_nexp}
\deqn{mean\_se\_up\_nexp = plot\_ancova\_mean\_se\_up\_nexp - plot\_ancova\_mean\_nexp}
\deqn{mean\_se\_nexp = \frac{mean\_se\_lo\_nexp + mean\_se\_up\_nexp}{2}}

Then, calculations of the \code{\link{es_from_ancova_means_se}} are used.
\enumerate{
\item If the bounds of the 95\% confidence intervals are provided, the calculations
of the \code{\link{es_from_ancova_means_ci}()} are used.
}
}
\examples{
es_from_plot_ancova_means(
  n_exp = 35, n_nexp = 35,
  cov_outcome_r = 0.2, n_cov_ancova = 4,
  plot_ancova_mean_exp = 89, plot_ancova_mean_nexp = 104,
  plot_ancova_mean_sd_lo_exp = 109, plot_ancova_mean_sd_lo_nexp = 83,
  plot_ancova_mean_sd_up_exp = 69, plot_ancova_mean_sd_up_nexp = 125
)
}
