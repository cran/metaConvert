#' Convert a 2x2 table into several effect size measures
#'
#' @param n_cases_exp number of cases/events in the exposed group
#' @param n_cases_nexp number of cases/events in the non exposed group
#' @param n_controls_exp number of controls/no-event in the exposed group
#' @param n_controls_nexp number of controls/no-event in the non exposed group
#' @param table_2x2_to_cor formula used to obtain a correlation coefficient from the contingency table (see details).
#' @param reverse_2x2 a logical value indicating whether the direction of generated effect sizes should be flipped.
#'
#' @details
#' This function first computes (log) odds ratio (OR), (log) risk ratio (RR) and number needed to treat (NNT)
#' from the 2x2 table.
#' Cohen's d (D), Hedges' g (G) and correlation coefficients (R/Z) are then estimated from the OR.
#'
#' **To estimate an OR**, the formulas used (Box 6.4.a in the Cochrane Handbook) are:
#' \deqn{logor = log(\frac{n\_cases\_exp / n\_cases\_nexp}{n\_controls\_exp / n\_controls\_nexp})}
#' \deqn{logor\_se = \sqrt{\frac{1}{n\_cases\_exp} + \frac{1}{n\_cases\_nexp} + \frac{1}{n\_controls\_exp} + \frac{1}{n\_controls\_nexp}}}
#'
#' **To estimate an RR**, the formulas used (Box 6.4.a in the Cochrane Handbook) are:
#' \deqn{logrr = log(\frac{n\_cases\_exp / n\_exp}{n\_cases\_nexp / n\_nexp})}
#' \deqn{logrr\_se = \sqrt{\frac{1}{n\_cases\_exp} - \frac{1}{n\_exp} + \frac{1}{n\_cases\_nexp} - \frac{1}{n\_nexp}}}
#'
#' **To estimate a NNT**, the formulas used are (Sedwick, 2013) :
#' \deqn{pt = \frac{n\_cases\_exp}{n\_cases\_exp + n\_controls\_exp}}
#' \deqn{pc = \frac{n\_cases\_nexp}{n\_cases\_nexp + n\_controls\_nexp}}
#' \deqn{AAR = pc - pt}
#' \deqn{nnt = \frac{1}{AAR}}
#'
#' **To convert the 2x2 table into a SMD**,
#' the function estimates an OR value from the 2x2 table (formula above)
#' that is then converted to a SMD
#' (see formula in \code{\link{es_from_or_se}()}).
#'
#' **To convert the 2x2 table into a correlation coefficient**,
#' various formulas can be used.
#'
#' **A.** First, Cooper et al. (2019) - \code{table_2x2_to_cor = "cooper"} -
#' proposes to convert the
#' 2x2 table into a OR (formula above), to convert this OR into a SMD
#' (see formula in \code{\link{es_from_or_se}()}), and to convert this
#' SMD into a correlation coefficient (see formula in \code{\link{es_from_cohen_d}()},
#' with the option \code{"smd_to_cor = 'lipsey_cooper'"}).
#'
#' **B.** Second, a correlation coefficient (more precisely - a phi coefficient)
#' can be obtained from the contingency table using the formula given in
#' Lipsey and Wilson (2001) - \code{table_2x2_to_cor = "lipsey"}.
#' The formulas used to estimate the r and z are:
#' \deqn{r = \frac{(n\_cases\_exp*n\_controls\_nexp - n\_controls\_exp*n\_cases\_nexp)}{\sqrt{(n\_exp) * (n\_nexp) * (n\_cases) * (n\_controls\_exp+n\_cases\_nexp)}}}
#' \deqn{z = atanh(r)}
#' \deqn{z\_se = logor\_se^2 * \frac{z^2}{\log(or)^2}}
#' \deqn{z\_ci\_lo = z - qnorm(.975)*z\_se}
#' \deqn{z\_ci\_up = z + qnorm(.975)*z\_se}
#' \deqn{r\_ci\_lo = tanh(z\_ci\_lo)}
#' \deqn{r\_ci\_up = tanh(z\_ci\_up)}
#' \deqn{effective\_n = \frac{1}{z\_se^2 + 3}}
#' \deqn{r\_se = \frac{(1 - r^2)^2}{effective\_n - 1}}
#'
#' **C.** Third, a tetrachoric correlation can be estimated from the 2x2 table
#' - \code{table_2x2_to_cor = "tetrachoric"}.
#' Given the heavy calculations required for this effect size measure,
#' we relied on the implementation of the formulas of the 'metafor' package. More
#' information can be retrieved here
#' (https://wviechtb.github.io/metafor/reference/escalc.html#-b-measures-for-two-dichotomous-variables).
#'
#' @return
#' This function estimates and converts between several effect size measures.
#'
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab OR + RR + NNT\cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab D + G + R + Z \cr
#'  \tab \cr
#'  \code{required input data} \tab See 'Section 7. Contingency (2x2) table or proportions'\cr
#'  \tab https://metaconvert.org/html/input.html\cr
#'  \tab \cr
#' }
#'
#' @md
#'
#' @export es_from_2x2
#'
#' @references
#' Cooper, H., Hedges, L.V., & Valentine, J.C. (Eds.). (2019). The handbook of research synthesis and meta-analysis. Russell Sage Foundation.
#'
#' Higgins JPT, Thomas J, Chandler J, Cumpston M, Li T, Page MJ, Welch VA (editors). Cochrane Handbook for Systematic Reviews of Interventions version 6.3 (updated February 2022). Available from www.training.cochrane.org/handbook.
#'
#' Lipsey, M. W., & Wilson, D. B. (2001). Practical meta-analysis. Sage Publications, Inc.
#'
#' Sedgwick, P. (2013). What is number needed to treat (NNT)? Bmj, 347.
#'
#' @examples
#' es_from_2x2(n_cases_exp = 467, n_cases_nexp = 22087, n_controls_exp = 261, n_controls_nexp = 8761)
es_from_2x2 <- function(n_cases_exp, n_cases_nexp,
                        n_controls_exp, n_controls_nexp,
                        table_2x2_to_cor = "tetrachoric", reverse_2x2) {
  if (missing(reverse_2x2)) reverse_2x2 <- rep(FALSE, length(n_cases_exp))
  reverse_2x2[is.na(reverse_2x2)] <- FALSE
  if (length(reverse_2x2) == 1) reverse_2x2 = c(rep(reverse_2x2, length(n_cases_exp)))
  if (length(reverse_2x2) != length(n_cases_exp)) stop("The length of the 'reverse_2x2' argument of incorrectly specified.")

  if (!all(table_2x2_to_cor %in% c("tetrachoric", "cooper", "lipsey"))) {
    stop(paste0(
      "'",
      unique(table_2x2_to_cor[!table_2x2_to_cor %in% c("tetrachoric", "cooper", "lipsey")]),
      "' not in tolerated values for the 'table_2x2_to_cor' argument.",
      " Possible inputs are: 'tetrachoric', 'cooper', 'lipsey'"
    ))
  }

  zero <- which(n_cases_exp == 0 | n_cases_nexp == 0 | n_controls_exp == 0 | n_controls_nexp == 0)
  n_cases_exp[zero] <- n_cases_exp[zero] + 0.5
  n_cases_nexp[zero] <- n_cases_nexp[zero] + 0.5
  n_controls_exp[zero] <- n_controls_exp[zero] + 0.5
  n_controls_nexp[zero] <- n_controls_nexp[zero] + 0.5
  n_exp <- n_cases_exp + n_controls_exp
  n_nexp <- n_cases_nexp + n_controls_nexp

  or_raw <- suppressWarnings((n_cases_exp * n_controls_nexp) / (n_controls_exp * n_cases_nexp))
  or <- ifelse(reverse_2x2, 1 / or_raw, or_raw)
  se_or <- suppressWarnings(sqrt(1 / n_cases_exp + 1 / n_controls_exp + 1 / n_cases_nexp + 1 / n_controls_nexp))

  rr_raw <- suppressWarnings((n_cases_exp / n_exp) / (n_cases_nexp / n_nexp))
  rr <- ifelse(reverse_2x2, 1 / rr_raw, rr_raw)
  se_rr <- suppressWarnings(sqrt(1 / n_cases_exp - 1 / n_exp + 1 / n_cases_nexp - 1 / n_nexp))

  logOR <- suppressWarnings(log(or))
  d <- logOR * sqrt(3) / pi
  d_se <- sqrt(se_or^2 * 3 / (pi^2))

  es <- .es_from_d(
    d = d, d_se = d_se,
    n_exp = n_exp,
    n_nexp = n_nexp,
    smd_to_cor = rep("lipsey_cooper", length(d))
  )

  es$logor <- logOR
  es$logor_se <- se_or
  es$logor_ci_lo <- es$logor - qnorm(.975) * es$logor_se
  es$logor_ci_up <- es$logor + qnorm(.975) * es$logor_se

  es$logrr <- log(rr)
  es$logrr_se <- se_rr
  es$logrr_ci_lo <- es$logrr - qnorm(.975) * es$logrr_se
  es$logrr_ci_up <- es$logrr + qnorm(.975) * es$logrr_se

  dat2x2 <- data.frame(
    n_cases_exp = n_cases_exp, n_controls_exp = n_controls_exp,
    n_cases_nexp = n_cases_nexp, n_controls_nexp = n_controls_nexp,
    n_exp = n_exp, n_nexp = n_nexp,
    table_2x2_to_cor = table_2x2_to_cor, reverse_2x2 = reverse_2x2
  )

  nn_miss <- which(
    (dat2x2$table_2x2_to_cor == "lipsey" |
       dat2x2$table_2x2_to_cor == "tetrachoric") &
      !is.na(dat2x2$n_cases_exp) & !is.na(dat2x2$n_controls_exp) &
      !is.na(dat2x2$n_cases_nexp) & !is.na(dat2x2$n_controls_nexp)
  )

  if (length(nn_miss) != 0) {
    no_cooper = which(dat2x2$table_2x2_to_cor != "cooper")

    es$r[no_cooper] <- es$r_se[no_cooper] <-
      es$r_ci_lo[no_cooper] <- es$r_ci_up[no_cooper] <-
      es$z[no_cooper] <- es$z_se[no_cooper] <-
      es$z_ci_lo[no_cooper] <- es$z_ci_up[no_cooper] <- NA

    res_tet <- suppressWarnings(t(mapply(.contigency_to_cor,
      n_cases_exp = dat2x2$n_cases_exp[nn_miss],
      n_controls_exp = dat2x2$n_controls_exp[nn_miss],
      n_cases_nexp = dat2x2$n_cases_nexp[nn_miss],
      n_controls_nexp = dat2x2$n_controls_nexp[nn_miss],
      table_2x2_to_cor = dat2x2$table_2x2_to_cor[nn_miss],
      reverse_2x2 = dat2x2$reverse_2x2[nn_miss]
    )))

    es$r[nn_miss] <- res_tet[, 1]
    es$r_se[nn_miss] <- suppressWarnings(sqrt(res_tet[, 2]))
    es$r_ci_lo[nn_miss] <- res_tet[, 3]
    es$r_ci_up[nn_miss] <- res_tet[, 4]

    es$z[nn_miss] <- res_tet[, 5]
    es$z_se[nn_miss] <- suppressWarnings(sqrt(res_tet[, 6]))
    es$z_ci_lo[nn_miss] <- res_tet[, 7]
    es$z_ci_up[nn_miss] <- res_tet[, 8]
  }

  pc <- n_cases_nexp / (n_cases_nexp+ n_controls_nexp)
  pt <- n_cases_exp / (n_cases_exp + n_controls_exp )
  AAR <- pc - pt
  es$nnt <- ifelse(reverse_2x2, -1/AAR, 1/AAR)

  es$info_used <- "2x2"
  return(es)
}

#' Convert a table with the number of cases and row marginal sums into several effect size measures
#'
#' @param n_cases_exp number of cases/events in the exposed group
#' @param n_cases_nexp number of cases/events in the non exposed group
#' @param n_exp total number of participants in the exposed group
#' @param n_nexp total number of participants in the non exposed group
#' @param table_2x2_to_cor formula used to obtain a correlation coefficient from the contigency table (see details).
#' @param reverse_2x2 a logical value indicating whether the direction of generated effect sizes should be flipped.
#'
#' @details
#' This function uses the number of cases in both the exposed
#' and non-exposed groups and the total number of participants exposed and non-exposed
#' to recreate a 2x2 table.
#' Then relies on the calculations of the \code{\link{es_from_2x2}} function.
#' \deqn{n\_controls\_exp = n\_exp - n\_cases\_exp}
#' \deqn{n\_controls\_nexp = n\_nexp - n\_cases\_nexp}
#'
#' @return
#' This function estimates and converts between several effect size measures.
#'
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab OR + RR + NNT\cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab D + G + R + Z \cr
#'  \tab \cr
#'  \code{required input data} \tab See 'Section 7. Contingency (2x2) table or proportions'\cr
#'  \tab https://metaconvert.org/html/input.html\cr
#'  \tab \cr
#' }
#'
#'
#' @md
#'
#' @export es_from_2x2_sum
#'
#' @examples
#' es_from_2x2_sum(n_cases_exp = 10, n_exp = 40, n_cases_nexp = 25, n_nexp = 47)
es_from_2x2_sum <- function(n_cases_exp, n_exp, n_cases_nexp, n_nexp,
                            table_2x2_to_cor = "tetrachoric", reverse_2x2) {
  if (missing(reverse_2x2)) reverse_2x2 <- rep(FALSE, length(n_cases_exp))
  reverse_2x2[is.na(reverse_2x2)] <- FALSE

  es <- es_from_2x2(
    n_cases_exp = n_cases_exp,
    n_cases_nexp = n_cases_nexp,
    n_controls_exp = n_exp - n_cases_exp,
    n_controls_nexp = n_nexp - n_cases_nexp,
    table_2x2_to_cor = table_2x2_to_cor,
    reverse_2x2 = reverse_2x2
  )

  es$info_used <- "2x2_sum"

  return(es)
}

#' Convert the proportion of occurrence of a binary event in two independent groups into several effect size measures
#'
#' @param prop_cases_exp proportion of cases/events in the exposed group (ranging from 0 to 1)
#' @param prop_cases_nexp proportion of cases/events in the non-exposed group (ranging from 0 to 1)
#' @param n_exp total number of participants in the exposed group
#' @param n_nexp total number of participants in the non exposed group
#' @param table_2x2_to_cor formula used to obtain a correlation coefficient from the contigency table (see details).
#' @param reverse_prop a logical value indicating whether the direction of generated effect sizes should be flipped.
#'
#'
#' @details
#' This function uses the proportions and sample size to
#' recreate the 2x2 table, and
#' then relies on the calculations of the \code{\link{es_from_2x2_sum}()} function.
#'
#' The formulas used is to obtain the 2x2 table are
#' \deqn{n\_cases\_exp = prop\_cases\_exp * n\_exp}
#' \deqn{n\_cases\_nexp = prop\_cases\_nexp * n\_nexp}
#' \deqn{n\_controls\_exp = (1 - prop\_cases\_exp) * n\_exp}
#' \deqn{n\_controls\_nexp = (1 - prop\_cases\_nexp) * n\_nexp}
#'
#' @return
#' This function estimates and converts between several effect size measures.
#'
#' \tabular{ll}{
#'  \code{natural effect size measure} \tab OR + RR + NNT\cr
#'  \tab \cr
#'  \code{converted effect size measure} \tab D + G + R + Z \cr
#'  \tab \cr
#'  \code{required input data} \tab See 'Section 7. Contingency (2x2) table or proportions'\cr
#'  \tab https://metaconvert.org/html/input.html\cr
#'  \tab \cr
#' }
#'
#' @export es_from_2x2_prop
#'
#'
#' @md
#'
#' @examples
#' es_from_2x2_prop(prop_cases_exp = 0.80, prop_cases_nexp = 0.60, n_exp = 10, n_nexp = 20)
es_from_2x2_prop <- function(prop_cases_exp, prop_cases_nexp, n_exp, n_nexp,
                             table_2x2_to_cor = "tetrachoric", reverse_prop) {
  if (missing(reverse_prop)) reverse_prop <- rep(FALSE, length(prop_cases_exp))
  reverse_prop[is.na(reverse_prop)] <- FALSE

  n_cases_exp <- round(prop_cases_exp * n_exp)
  n_cases_nexp <- round(prop_cases_nexp * n_nexp)
  n_controls_exp <- n_exp - n_cases_exp
  n_controls_nexp <- n_nexp - n_cases_nexp

  es <- es_from_2x2(
    n_cases_exp = n_cases_exp, n_cases_nexp = n_cases_nexp,
    n_controls_exp = n_controls_exp, n_controls_nexp = n_controls_nexp,
    table_2x2_to_cor = table_2x2_to_cor, reverse_2x2 = reverse_prop
  )

  es$info_used <- "2x2_prop"

  return(es)
}

