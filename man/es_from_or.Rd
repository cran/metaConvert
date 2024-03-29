% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/es_from_stand_OR.R
\name{es_from_or}
\alias{es_from_or}
\title{Convert an odds ratio value to several effect size measures}
\usage{
es_from_or(
  or,
  logor,
  n_cases,
  n_controls,
  n_sample,
  small_margin_prop,
  baseline_risk,
  n_exp,
  n_nexp,
  or_to_cor = "bonett",
  or_to_rr = "metaumbrella_cases",
  reverse_or
)
}
\arguments{
\item{or}{odds ratio value}

\item{logor}{log odds ratio value}

\item{n_cases}{number of cases/events}

\item{n_controls}{number of controls/no-event}

\item{n_sample}{total number of participants in the sample}

\item{small_margin_prop}{smallest margin proportion of the underlying 2x2 table}

\item{baseline_risk}{proportion of cases in the non-exposed group}

\item{n_exp}{number of participants in the exposed group}

\item{n_nexp}{number of participants in the non-exposed group}

\item{or_to_cor}{formula used to convert the \code{or} value into a correlation coefficient (see details).}

\item{or_to_rr}{formula used to convert the \code{or} value into a risk ratio (see details).}

\item{reverse_or}{a logical value indicating whether the direction of the generated effect sizes should be flipped.}
}
\value{
This function estimates and converts between several effect size measures.

\tabular{ll}{
\code{natural effect size measure} \tab N/A\cr
\tab \cr
\code{converted effect size measure} \tab OR + RR + NNT\cr
\code{} \tab D + G + R + Z\cr
\tab \cr
\code{required input data} \tab See 'Section 2. Odds Ratio'\cr
\tab https://metaconvert.org/html/input.html\cr
\tab \cr
}
}
\description{
Convert an odds ratio value to several effect size measures
}
\details{
This function computes the standard error of the log odds ratio.
Risk ratio (RR), Cohen's d (D), Hedges' g (G) and correlation coefficients (R/Z),
are converted from the odds ratio value.

\strong{Estimation of the standard error of the log OR.}
This function generates the standard error of an odds ratio (OR) based on the OR value and the number of cases and controls.
More precisely, this function simulates all combinations of the possible number of cases and controls in the exposed and non-exposed groups
compatible with the reported OR value and with the overall number of cases and controls. Then, our function assumes that the variance of the OR
is equal to the mean of the standard error of all possible situations. This estimation thus necessarily comes with some imprecision and should not
be used before having requested the value (or raw data) to authors of the original report.

\strong{Conversion of other effect size measures.}
Calculations of \code{\link{es_from_or_se}()} are then applied to
estimate the other effect size measures
}
\examples{
es_or_guess <- es_from_or(or = 0.5, n_cases = 210, n_controls = 220)
es_or <- es_from_or_se(or = 0.5, logor_se = 0.4, n_cases = 210, n_controls = 220)
round(es_or_guess$logor_se, 0.10) == round(es_or$logor_se, 0.10)
}
\references{
Gosling, C. J., Solanes, A., Fusar-Poli, P., & Radua, J. (2023). metaumbrella: the first comprehensive suite to perform data analysis in umbrella reviews with stratification of the evidence. BMJ mental health, 26(1), e300534. https://doi.org/10.1136/bmjment-2022-300534
}
