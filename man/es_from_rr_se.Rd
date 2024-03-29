% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/es_from_stand_RR.R
\name{es_from_rr_se}
\alias{es_from_rr_se}
\title{Convert a risk ratio value and standard error to three effect measures (SMD, OR, COR)}
\usage{
es_from_rr_se(
  rr,
  logrr,
  logrr_se,
  baseline_risk,
  n_exp,
  n_nexp,
  n_cases,
  n_controls,
  smd_to_cor = "viechtbauer",
  rr_to_or = "metaumbrella",
  reverse_rr
)
}
\arguments{
\item{rr}{risk ratio value}

\item{logrr}{log risk ratio value}

\item{logrr_se}{standard error of the log risk ratio}

\item{baseline_risk}{proportion of cases in the non-exposed group}

\item{n_exp}{number of participants in the exposed group}

\item{n_nexp}{number of participants in the non-exposed group}

\item{n_cases}{number of cases/events}

\item{n_controls}{number of controls/no-event}

\item{smd_to_cor}{formula used to convert the SMD value (converted from RR) into a coefficient correlation (see \code{\link{es_from_cohen_d}}).}

\item{rr_to_or}{formula used to convert the \code{rr} value into an odds ratio (see details).}

\item{reverse_rr}{a logical value indicating whether the direction of the generated effect sizes should be flipped.}
}
\value{
This function estimates and converts between several effect size measures.

\tabular{ll}{
\code{natural effect size measure} \tab RR\cr
\tab \cr
\code{converted effect size measure} \tab OR + NNT\cr
\tab \cr
\code{required input data} \tab See 'Section 3. Risk Ratio'\cr
\tab https://metaconvert.org/html/input.html\cr
\tab \cr
}
}
\description{
Convert a risk ratio value and standard error to three effect measures (SMD, OR, COR)
}
\details{
This function converts the (log) risk ratio (RR) value and its standard error
to odds ratio (OR) and number needed to treat.

\strong{To estimate the odds ratio and its standard error}, various formulas can be used.

\strong{A.} First, the approach described in Grant (2014) can be used.
However, in the paper, only the formula to convert an RR value to a OR value
is described. To derive the variance, we used this formula to convert the bounds of the 95\% CI, which
were then used to obtain the variance.

This argument requires (rr + baseline_risk + rr_ci_lo + rr_ci_up) to generate a RR.
The following formulas are used (br = baseline_risk):
\deqn{or = \frac{rr * (1 - br)}{1 - rr * br}}
\deqn{or\_ci\_lo = \frac{rr\_ci\_lo}{1 - br + br*rr\_ci\_lo}}
\deqn{or\_ci\_up = \frac{rr\_ci\_up}{1 - br + br*rr\_ci\_up}}
\deqn{logor\_se = \frac{log(or\_ci\_up) - log(or\_ci\_lo)}{2 * qnorm(.975)}}

\strong{B.} Second, the formulas implemented in the metaumbrella package can be used (\code{or_to_rr = "metaumbrella_exp"}).
This argument requires (rr + logrr_se + n_exp + n_nexp) to generate a OR.
More precisely, we previously developed functions that simulate all combinations of the possible number of cases and controls
in the exposed and non-exposed groups compatible with the actual value of the RR.
Then, the functions select the contingency table whose standard error coincides best with the standard error reported.
The RR value and its standard are obtained from this estimated contingency table.

\strong{C.} Third, it is possible to transpose the RR to a OR (\code{rr_to_or = "transpose"}).
This argument requires (rr + logrr_se) to generate a OR.
It is known that OR and RR are similar when the baseline risk is small.
Therefore, users can request to simply transpose the RR value & standard error into a OR value & standard error.
\deqn{or = rr}
\deqn{logor\_se = logrr\_se}

\strong{D.} Fourth, it is possible to recreate the 2x2 table using the dipietrantonj's formulas (\code{rr_to_or = "dipietrantonj"}).
This argument requires (rr + logrr_ci_lo + logrr_ci_lo) to generate a OR. Information on this approach can be retrieved in
Di Pietrantonj (2006).

\strong{To estimate the NNT}, the formulas used are :
\deqn{nnt = \frac{1}{br * (1 - rr)}}

\strong{To estimate the Cohen's d value and its standard error}, the function first converts the RR value and standard error into OR and standard error,
and then converts these values into Cohen's d using the following formulas:
\deqn{cohen\_d = \log(or) * \frac{\sqrt{3}}{\pi}}
\deqn{cohen\_d\_se = \sqrt{\frac{logor\_se^2 * 3}{\pi^2}}}
}
\examples{
es_from_rr_se(rr = 2.12, logrr_se = 0.242, n_exp = 120, n_nexp = 44)
}
\references{
Di Pietrantonj C. (2006). Four-fold table cell frequencies imputation in meta analysis. Statistics in medicine, 25(13), 2299–2322. https://doi.org/10.1002/sim.2287

Gosling, C. J., Solanes, A., Fusar-Poli, P., & Radua, J. (2023). metaumbrella: the first comprehensive suite to perform data analysis in umbrella reviews with stratification of the evidence. BMJ mental health, 26(1), e300534. https://doi.org/10.1136/bmjment-2022-300534

Grant R. L. (2014). Converting an odds ratio to a range of plausible relative risks for better communication of research findings. BMJ (Clinical research ed.), 348, f7450. https://doi.org/10.1136/bmj.f7450

Veroniki, A. A., Pavlides, M., Patsopoulos, N. A., & Salanti, G. (2013). Reconstructing 2x2 contingency tables from odds ratios using the Di Pietrantonj method: difficulties, constraints and impact in meta-analysis results. Research synthesis methods, 4(1), 78–94. https://doi.org/10.1002/jrsm.1061
}
