#' metaConvert: An R Package Dedicated to Automated Effect Size Calculations
#'
#' The \pkg{metaConvert} package automatically estimates 11 effect size measures from a well-formatted dataframe.
#' Various other functions can help, for example, removing dependency between several effect sizes,
#' or identifying differences between two dataframes.
#' This package is mainly designed to assist in conducting a systematic review with a meta-analysis, but it can be
#' useful to any researcher interested in estimating an effect size.
#'
#' # Overview of the package
#' To visualize all the types of input data that can be used to estimate the 11 effect size measures available
#' in metaConvert, you can use the \code{\link{see_input_data}()} function.
#'
#' # Estimate effect sizes
#' To automatically estimate effect sizes directly from a dataset, you can use the \code{\link{convert_df}()} function.
#'
#' # Aggregate dependent effect sizes
#' To automatically aggregate dependent effect sizes using Borenstein's formulas,
#' you can use the \code{\link{aggregate_df}()} function. This function can handle dependent effect sizes
#' from multiple subgroups, or dependent effect sizes from the same participants.
#'
#' # Flag differences between two datasets
#' If pairs of data extractors have generated similar datasets that should be compared,
#' you can use the \code{\link{compare_df}()} function.
#'
#' # Prepare a dataset extraction sheet
#' If you have not started data extraction yet,
#' you can use the \code{\link{data_extraction_sheet}()} function to obtain a
#' perfectly formatted data extraction sheet.
#'
#' @details
#'
#' # Well-formatted dataset
#' One of the specificities of the \pkg{metaConvert} package is that its core function (\code{\link{convert_df}})
#' does not have arguments to specify the names of the variables contained in the dataset.
#' While this allow using a convenient automatic process in the calculations, this requires that the datasets
#' passed to this function respect a very precise formatting (which we will refer to as \code{well-formatted dataset}).
#'
#' Rather than a long description of all column names, we built several tools
#' that help you find required information.
#' 1. You can use the \code{\link{data_extraction_sheet}()} function that
#' generates an excel/csv/txt file containing all the column names available,
#' as well as a description of the information it should contain.
#' 2. You can use the \code{\link{see_input_data}()} function that generates
#' a list of all available types of input data as well as their estimated/converted
#' effect size measures. This function also points out to the corresponding helper tables
#' available in https://metaconvert.org
#'
#'
#' # Effect size measures available
#'
#' Eleven effect size measures are accepted:
#'
#' \itemize{
#'  \item \bold{"d"}: standardized mean difference (i.e., Cohen's d)
#'  \item \bold{"g"}: Hedges' g
#'  \item \bold{"md"}: mean difference
#'  \item \bold{"r"}: Correlation coefficient
#'  \item \bold{"z"}: Fisher's r-to-z correlation
#'  \item \bold{"or"} or \bold{"logor"}: odds ratio  or its logarithm
#'  \item \bold{"rr"} or \bold{"logrr"}: risk ratio or its logarithm
#'  \item \bold{"irr"} or \bold{"logirr"}: incidence rate ratio or its logarithm
#'  \item \bold{"nnt"}: number needed to treat
#'  \item \bold{"logcvr"}: log coefficient of variation
#'  \item \bold{"logvr"}: log variability ratio
#' }
#'
#' # Output
#' All the functions of the \pkg{metaConvert} package that are dedicated to effect size calculations
#' (i.e., all the functions named \code{es_from_*}) return a dataframe that contain,
#' depending on the function - some of the following columns:
#' \tabular{ll}{
#'  \code{info_used} \tab input data used to generate the effect size.\cr
#'  \tab \cr
#'  \code{md} \tab value of the mean difference.\cr
#'  \tab \cr
#'  \code{md_se} \tab standard error of the mean difference.\cr
#'  \tab \cr
#'  \code{md_ci_lo} \tab lower bound of the 95% CI of the mean difference.\cr
#'  \tab \cr
#'  \code{md_ci_up} \tab upper bound of the 95% CI of the mean difference.\cr
#'  \tab \cr
#'  \code{d} \tab value of the Cohen's d.\cr
#'  \tab \cr
#'  \code{d_se} \tab standard error of the Cohen's d.\cr
#'  \tab \cr
#'  \code{d_ci_lo} \tab lower bound of the 95% CI of the Cohen's d.\cr
#'  \tab \cr
#'  \code{d_ci_up} \tab upper bound of the 95% CI of the Cohen's d.\cr
#'  \tab \cr
#'  \code{g} \tab value of the Hedges' g.\cr
#'  \tab \cr
#'  \code{g_se} \tab standard error of the Hedges' g.\cr
#'  \tab \cr
#'  \code{g_ci_lo} \tab lower bound of the 95% CI of the Hedges' g.\cr
#'  \tab \cr
#'  \code{g_ci_up} \tab upper bound of the 95% CI of the Hedges' g.\cr
#'  \tab \cr
#'  \code{r} \tab value of the correlation coefficient.\cr
#'  \tab \cr
#'  \code{r_se} \tab standard error of the correlation coefficient.\cr
#'  \tab \cr
#'  \code{r_ci_lo} \tab lower bound of the 95% CI of the correlation coefficient.\cr
#'  \tab \cr
#'  \code{r_ci_up} \tab upper bound of the 95% CI of the correlation coefficient.\cr
#'  \tab \cr
#'  \code{z} \tab value of the r-to-z transformed correlation coefficient.\cr
#'  \tab \cr
#'  \code{z_se} \tab standard error of the r-to-z transformed correlation coefficient.\cr
#'  \tab \cr
#'  \code{z_ci_lo} \tab lower bound of the 95% CI of the r-to-z transformed correlation coefficient.\cr
#'  \tab \cr
#'  \code{z_ci_up} \tab upper bound of the 95% CI of the r-to-z transformed correlation coefficient.\cr
#'  \tab \cr
#'  \code{logor} \tab value of the log odds ratio.\cr
#'  \tab \cr
#'  \code{logor_se} \tab standard error of the log odds ratio.\cr
#'  \tab \cr
#'  \code{logor_ci_lo} \tab lower bound of the 95% CI of the log odds ratio.\cr
#'  \tab \cr
#'  \code{logor_ci_up} \tab upper bound of the 95% CI of the log odds ratio.\cr
#'  \tab \cr
#'  \code{logrr} \tab value of the log risk ratio.\cr
#'  \tab \cr
#'  \code{logrr_se} \tab standard error of the log risk ratio.\cr
#'  \tab \cr
#'  \code{logrr_ci_lo} \tab lower bound of the 95% CI of the log risk ratio.\cr
#'  \tab \cr
#'  \code{logrr_ci_up} \tab upper bound of the 95% CI of the log risk ratio.\cr
#'  \tab \cr
#'  \code{logirr} \tab value of the log incidence rate ratio.\cr
#'  \tab \cr
#'  \code{logirr_se} \tab standard error of the log incidence rate ratio.\cr
#'  \tab \cr
#'  \code{logirr_ci_lo} \tab lower bound of the 95% CI of the log incidence rate ratio.\cr
#'  \tab \cr
#'  \code{logirr_ci_up} \tab upper bound of the 95% CI of the log incidence rate ratio.\cr
#'  \tab \cr
#'  \code{logvr} \tab value of the log variability ratio.\cr
#'  \tab \cr
#'  \code{logvr_se} \tab standard error of the log variability ratio.\cr
#'  \tab \cr
#'  \code{logvr_ci_lo} \tab lower bound of the 95% CI of the log variability ratio.\cr
#'  \tab \cr
#'  \code{logvr_ci_up} \tab upper bound of the 95% CI of the log variability ratio.\cr
#'  \tab \cr
#'  \code{logcvr} \tab value of the log coefficient of variation.\cr
#'  \tab \cr
#'  \code{logcvr_se} \tab standard error of the log coefficient of variation.\cr
#'  \tab \cr
#'  \code{logcvr_ci_lo} \tab lower bound of the 95% CI of the log coefficient of variation.\cr
#'  \tab \cr
#'  \code{logcvr_ci_up} \tab upper bound of the 95% CI of the log coefficient of variation.\cr
#'  \tab \cr
#'  \code{nnt} \tab number needed to treat.\cr
#'  \tab \cr
#' }
#'
#' @md
#'
#' @docType package
#'
#' @name metaConvert-package
NULL

#' @importFrom stats sd aggregate.data.frame df dnorm na.omit optimize pnorm pt qnorm qt
NULL

utils::globalVariables(c(
  "y", "discard", "reshape",
  "row_id", "study_id", "author", "year", "predictor",
  "outcome", "info_expected", "adjusted_input", "df.haza"
))
