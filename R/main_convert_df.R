#' Automatically compute effect sizes from a well formatted dataset
#'
#' @param x a well formatted dataset
#' @param measure the effect size measure that will be estimated from the information stored in the dataset. See details.
#' @param main_es a logical variable indicating whether a main effect size should be selected when overlapping data are present. See details.
#' @param split_adjusted a logical value indicating whether crude and adjusted effect sizes should be presented separately. See details.
#' @param format_adjusted presentation format of the adjusted effect sizes. See details.
#' @param hierarchy a character string indicating the hierarchy in the information to be prioritized for the effect size calculations (see details). See details.
#' @param es_selected the method used to select the main effect size when several information allows to estimate an effect size for the same association/comparison. Must be either "minimum" (the smallest effect size will be selected), "maximum" (the largest effect size will be selected) or "hierarchy" (the effect size computed from the information specified highest in the hierarchy will be selected). See details.
#' @param rr_to_or formula used to convert the \code{rr} value into an odds ratio.
#' @param or_to_rr formula used to convert the \code{or} value into a risk ratio.
#' @param or_to_cor formula used to convert the \code{or} value into a correlation coefficient.
#' @param table_2x2_to_cor formula used to obtain a correlation coefficient from the contingency table.
#' @param pre_post_to_smd formula used to obtain a SMD from pre/post means and SD of two independent groups.
#' @param r_pre_post pre-post correlation across the two groups (use this argument only if the precise correlation in each group is unknown)
#' @param smd_to_cor formula used to convert the \code{cohen_d} value into a coefficient correlation.
#' @param cor_to_smd formula used to convert a correlation coefficient value into a SMD.
#' @param yates_chisq a logical value indicating whether the Chi square has been performed using Yate's correction for continuity.
#' @param unit_type the type of unit for the \code{unit_increase_iv} argument. Must be either "sd" or "value" (see \code{\link{es_from_pearson_r}}).
#' @param verbose a logical variable indicating whether text outputs and messages should be generated. We recommend turning this option to FALSE only after having carefully read all the generated messages.
#'
#' @details
#' This function automatically computes or converts between 11 effect sizes
#' measures from any relevant type of input data stored in the
#' dataset you pass to this function.
#'
#' ## Effect size measures
#' Possible effect size measures are:
#' 1. Cohen's d ("d")
#' 2. Hedges' g ("g")
#' 3. mean difference ("md")
#' 4. (log) odds ratio ("or" and "logor")
#' 5. (log) risk ratio ("rr" and "logrr")
#' 6. (log) incidence rate ratio ("irr" and "logirr")
#' 7. correlation coefficient ("r")
#' 8. transformed r-to-z correlation coefficient ("z")
#' 9. log variability ratio ("logvr")
#' 10. log coefficient of variation ("logcvr")
#' 11. number needed to treat ("nnt")
#'
#' ## Computation of a main effect size
#' If you enter multiple types of input data
#' (e.g., means/sd of two groups and a student t-test value)
#' for the same comparison i.e., for the same row of the dataset,
#' the \code{convert_df()} function can have two behaviours.
#' If you set:
#' - \code{main_es = FALSE} the function will estimate all possible effect sizes from all
#' types of input data (which implies that if a comparison has **several types of input data**,
#' it will result in **multiple rows** in the dataframe returned by the function)
#' - \code{main_es = TRUE} the function will select one effect size per comparison
#' (which implies that if a comparison has **several types of input data**,
#' it will result in a **unique row** in the dataframe returned by the function)
#'
#' ## Selection of input data for the computation of the main effect size
#' If you choose to estimate one main effect size (i.e., by setting \code{main_es = TRUE}),
#' you have several options to select this main effect size.
#' If you set:
#' - \code{es_selected = "hierarchy"}: the main effect size will be selected by prioritizing
#' specific types of input data over other (see next section "Hierarchy").
#' - \code{es_selected = "minimum"}: the main effect size will be selected by selecting
#' the lowest effect size available.
#' - \code{es_selected = "maximum"}: the main effect size will be selected by selecting
#' the highest effect size available.
#'
#' ## Hierarchy
#' More than 70 different combinations of input data can be used to estimate an effect size.
#' You can retrieve the effect size measures estimated by each combination of input data
#' in the \code{\link{see_input_data}()} function.
#'
#' In the hierarchy argument, each type of input data should be separated by the symbol ">".
#' You can prioritize some types of input data by placing them at the begining of the
#' hierarchy argument.
#' For example, if you set:
#' - \code{hierarchy = "means_sd > means_se > student_t"}, the convert_df function will prioritize
#' the means + SD, then the means + SE, then the Student's t-test to estimate the main effect
#' size.
#' - \code{hierarchy = "2x2 > or_se > phi"}, the convert_df function will prioritize
#' the contigency table, then the odds ratio value + SE, then the phi coefficient to estimate
#' the main effect size.
#'
#' Importantly, if none of the types of input data indicated in the \code{hierarchy} argument
#' can be used to estimate the target effect size measure,
#' the \code{convert_df()} function will automatically try to use other types of input
#' data to estimate an effect size.
#'
#' ## Adjusted effect sizes
#' Some datasets will be composed of crude (i.e., non-adjusted) types of input data
#' (such as standard means + SD, Student's t-test, etc.) and adjusted types of input data
#' (such as means + SE from an ANCOVA model, a t-test from an ANCOVA, etc.).
#'
#' In these situations, you can decide to:
#' - treat crude and adjusted input data the same way \code{split_adjusted = FALSE}
#' - split calculations for crude and adjusted types of input data \code{split_adjusted = TRUE}
#'
#' If you want to split the calculations, you can decide to present the final dataset:
#' - in a long format (i.e., crude and adjusted effect sizes presented in separate rows \code{format_adjusted = "long"})
#' - in a wide format (i.e., crude and adjusted effect sizes presented in separate columns \code{format_adjusted = "wide"})
#'
#' @return
#' The \code{convert_df()} function returns a list of
#' more than 70 dataframes
#' (one for each function automatically applied to the dataset).
#' These dataframes systematically contain the columns described in
#' \code{\link{metaConvert-package}}.
#' The list of dataframes can be easily converted to a single,
#' calculations-ready dataframe
#' using the summary function (see \code{\link{summary.metaConvert}}).
#'
#' @md
#'
#' @export convert_df
#'
#' @examples
#' res <- convert_df(df.haza,
#'   measure = "g",
#'   split_adjusted = TRUE,
#'   es_selected = "minimum",
#'   format_adjusted = "long"
#' )
#' summary(res)
convert_df <- function(x, measure = c("d", "g", "md", "logor", "logrr", "logirr",
                                      "nnt", "r", "z", "logvr", "logcvr"),
                       main_es = TRUE,
                       es_selected = c("hierarchy", "minimum", "maximum"),
                       split_adjusted = TRUE,
                       format_adjusted = c("wide", "long"),
                       verbose = TRUE,
                       hierarchy = "means_sd > means_se > means_ci",
                       rr_to_or = "metaumbrella",
                       or_to_rr = "metaumbrella_cases",
                       or_to_cor = "bonett",
                       table_2x2_to_cor = "lipsey",
                       smd_to_cor = "viechtbauer",
                       pre_post_to_smd = "bonett",
                       r_pre_post = 0.5,
                       cor_to_smd = "viechtbauer",
                       unit_type = "raw_scale",
                       yates_chisq = FALSE) {

  # x = dat
  # measure = "d";
  # main_es = TRUE;
  # es_selected = "hierarchy";
  # split_adjusted = TRUE;
  # format_adjusted = "wide";
  # verbose = TRUE;
  # hierarchy = "means_sd > means_se > means_ci";
  # rr_to_or = "metaumbrella";
  # or_to_rr = "dipietrantonj";
  # or_to_cor = "bonett";
  # table_2x2_to_cor = "lipsey";
  # smd_to_cor = "viechtbauer";
  # chisq_to_cor = "tetrachoric";
  # phi_to_cor = "tetrachoric";
  # pre_post_to_smd = "bonett";
  # cor_to_smd = "viechtbauer";
  # unit_type = "raw_scale";
  # yates_chisq = FALSE
  measure = measure[1]
  es_selected = es_selected[1]
  format_adjusted = format_adjusted[1]

  x <- .check_data(x, split_adjusted = split_adjusted,
                   main_es = main_es, format = format_adjusted)
  x[is.na(x[, "r_pre_post_exp"]), "r_pre_post_exp"] <- r_pre_post
  x[is.na(x[, "r_pre_post_nexp"]), "r_pre_post_nexp"] <- r_pre_post
  r_pre_post = rep(r_pre_post, nrow(x))
  for (i in c("rr_to_or",
              "or_to_rr",
              "or_to_cor",
              "table_2x2_to_cor",
              "smd_to_cor",
              "pre_post_to_smd",
              "cor_to_smd",
              "unit_type")) {
    if (i %in% colnames(x)) {
      if (length(sapply(mget(i), function(x) x)) > 1) {
        stop("The length of the argument '", i, "' should be 1, or should be specified in a column of the dataset.")
      }
      x[, i][is.na(x[, i])] <- sapply(mget(i), function(x) x)
    } else {
      if (length(sapply(mget(i), function(x) x)) > 1) {
        stop("The length of the argument '", i, "' should be 1, or should be specified in a column of the dataset.")
      }
      x[, i] <- sapply(mget(i), function(x) x)
    }
  }

  if (verbose) message("Calculations in progress, it may take up to 30 sec...")
  measure = tolower(measure)
  if (!measure %in% c("d", "g", "md", "r", "z", "or", "rr", "irr", "logor", "logrr", "logirr", "logvr", "logcvr", "nnt")) {
    stop(paste0("'", measure, "' not in tolerated measures. Possible inputs are: 'md', 'd', 'g', 'or', 'rr', 'irr', 'logor', 'logrr', 'logirr', 'r', 'z', 'logvr', 'logcvr', 'nnt'"))
  } else if (!split_adjusted %in% c(TRUE, FALSE)) {
    stop(paste0("'", split_adjusted, "' not in tolerated values for the 'split_adjusted' argument. Should be a logical value (TRUE/FALSE)"))
  } else if (!format_adjusted %in% c("wide", "long")) {
    stop(paste0("'", format_adjusted, "' not in tolerated values for the 'format_adjusted' argument. Should be either 'long' or 'wide'."))
  } else if (!es_selected %in% c("minimum", "maximum", "hierarchy")) {
    stop(paste0("'", es_selected, "' not in tolerated values for the 'es_selected' argument. Possible inputs are: 'minimum', 'maximum', 'hierarchy'"))
  } else if (!main_es %in% c(TRUE, FALSE)) {
    stop(paste0("'", main_es, "' not in tolerated values for the 'main_es' argument. Should be a logical value (TRUE/FALSE)"))
  }

  if (measure %in% c("or", "rr", "irr")) {
    exp <- TRUE
    if (measure == "or") {
      measure <- "logor"
    } else if (measure == "rr") {
      measure <- "logrr"
    } else if (measure == "irr") {
      measure <- "logirr"
    }
  } else {
    exp <- FALSE
  }

  # SMD ----------------------------------------------
  es_cohen_d <- with(x, es_from_cohen_d(
    cohen_d = cohen_d, n_exp = n_exp, n_nexp = n_nexp,
    smd_to_cor = smd_to_cor, reverse_d = reverse_d)
  )

  es_cohen_d_adj <- with(x, es_from_cohen_d_adj(
    cohen_d_adj = cohen_d_adj, n_cov_ancova = n_cov_ancova, cov_outcome_r = cov_outcome_r,
    n_exp = n_exp, n_nexp = n_nexp, smd_to_cor = smd_to_cor, reverse_d = reverse_d
  ))
  es_hedges_g <- with(x, es_from_hedges_g(
    hedges_g = hedges_g, n_exp = n_exp, n_nexp = n_nexp,
    smd_to_cor = smd_to_cor, reverse_g = reverse_g
  ))


  # OR ----------------------------------------------
  es_odds_ratio <- with(x, es_from_or(
    or = or, logor = logor, baseline_risk = baseline_risk, n_exp = n_exp, n_nexp = n_nexp,
    n_cases = n_cases, n_controls = n_controls, n_sample = n_sample, small_margin_prop = small_margin_prop,
    or_to_rr = or_to_rr, or_to_cor = or_to_cor, reverse_or = reverse_or
  ))
  es_odds_ratio_se <- with(x, es_from_or_se(
    or = or, logor = logor, logor_se = logor_se, small_margin_prop = small_margin_prop,
    baseline_risk = baseline_risk,
    n_exp = n_exp, n_nexp = n_nexp, n_sample = n_sample,
    n_cases = n_cases, n_controls = n_controls,
    or_to_rr = or_to_rr, or_to_cor = or_to_cor,
    reverse_or = reverse_or
  ))
  es_odds_ratio_ci <- with(x, es_from_or_ci(
    or = or, or_ci_lo, or_ci_up, logor = logor, small_margin_prop = small_margin_prop,
    baseline_risk = baseline_risk, n_exp = n_exp, n_nexp = n_nexp,
    logor_ci_lo = logor_ci_lo, logor_ci_up = logor_ci_up,
    n_cases = n_cases, n_controls = n_controls, n_sample = n_sample,
    or_to_rr = or_to_rr, or_to_cor = or_to_cor, reverse_or = reverse_or
  ))
  es_odds_ratio_pval <- with(x, es_from_or_pval(
    or = or, logor = logor, or_pval = or_pval, small_margin_prop = small_margin_prop,
    baseline_risk = baseline_risk, n_exp = n_exp, n_nexp = n_nexp,
    n_cases = n_cases, n_controls = n_controls, n_sample = n_sample,
    or_to_rr = or_to_rr, or_to_cor = or_to_cor, reverse_or_pval = reverse_or_pval
  ))
  # R ----------------------------------------------
  es_pearson_r <- with(x, es_from_pearson_r(
    pearson_r = pearson_r, n_sample = n_sample,
    unit_type = unit_type, sd_iv = sd_iv,
    unit_increase_iv = unit_increase_iv,
    n_exp = n_exp, n_nexp = n_nexp, cor_to_smd = cor_to_smd,
    reverse_pearson_r = reverse_pearson_r
  ))
  es_fisher_z <- with(x, es_from_fisher_z(
    fisher_z = fisher_z, n_sample = n_sample,
    unit_type = unit_type, sd_iv = sd_iv,
    unit_increase_iv = unit_increase_iv,
    n_exp = n_exp, n_nexp = n_nexp, cor_to_smd = cor_to_smd,
    reverse_fisher_z = reverse_fisher_z
  ))

  # MEANS ----------------------------------------------
  es_means_sd_raw <- with(x, es_from_means_sd(
    n_exp = n_exp, n_nexp = n_nexp,
    mean_exp = mean_exp, mean_sd_exp = mean_sd_exp,
    mean_nexp = mean_nexp, mean_sd_nexp = mean_sd_nexp,
    smd_to_cor = smd_to_cor, reverse_means = reverse_means
  )
  )
  es_means_se_raw <- with(
    x,
    es_from_means_se(
      n_exp = n_exp, n_nexp = n_nexp,
      mean_exp = mean_exp, mean_se_exp = mean_se_exp,
      mean_nexp = mean_nexp, mean_se_nexp = mean_se_nexp,
      smd_to_cor = smd_to_cor, reverse_means = reverse_means
    )
  )
  es_means_ci_raw <- with(
    x,
    es_from_means_ci(
      n_exp = n_exp, n_nexp = n_nexp,
      mean_exp = mean_exp, mean_ci_lo_exp = mean_ci_lo_exp, mean_ci_up_exp = mean_ci_up_exp,
      mean_nexp = mean_nexp, mean_ci_lo_nexp = mean_ci_lo_nexp, mean_ci_up_nexp = mean_ci_up_nexp,
      smd_to_cor = smd_to_cor, reverse_means = reverse_means
    )
  )
  es_means_sd_pooled <- with(
    x,
    es_from_means_sd_pooled(
      n_exp = n_exp, n_nexp = n_nexp, mean_exp = mean_exp, mean_nexp = mean_nexp,
      mean_sd_pooled = mean_sd_pooled, smd_to_cor = smd_to_cor, reverse_means = reverse_means
    )
  )
  # plot
  es_plot_means_raw <- with(x, es_from_plot_means(
    n_exp = n_exp, n_nexp = n_nexp,
    plot_mean_exp = plot_mean_exp, plot_mean_nexp = plot_mean_nexp,
    plot_mean_sd_lo_exp = plot_mean_sd_lo_exp, plot_mean_sd_lo_nexp = plot_mean_sd_lo_nexp,
    plot_mean_sd_up_exp = plot_mean_sd_up_exp, plot_mean_sd_up_nexp = plot_mean_sd_up_nexp,
    plot_mean_se_lo_exp = plot_mean_se_lo_exp, plot_mean_se_lo_nexp = plot_mean_se_lo_nexp,
    plot_mean_se_up_exp = plot_mean_se_up_exp, plot_mean_se_up_nexp = plot_mean_se_up_nexp,
    plot_mean_ci_lo_exp = plot_mean_ci_lo_exp, plot_mean_ci_lo_nexp = plot_mean_ci_lo_nexp,
    plot_mean_ci_up_exp = plot_mean_ci_up_exp, plot_mean_ci_up_nexp = plot_mean_ci_up_nexp,
    smd_to_cor = smd_to_cor,
    reverse_plot_means = reverse_plot_means
  ))
  # PRE POST MEANS ----------------------------------------------
  es_means_sd_pre_post <- with(
    x,
    es_from_means_sd_pre_post(
      n_exp = n_exp, n_nexp = n_nexp,
      mean_pre_exp = mean_pre_exp, mean_exp = mean_exp,
      mean_pre_sd_exp = mean_pre_sd_exp, mean_sd_exp = mean_sd_exp,
      mean_pre_nexp = mean_pre_nexp, mean_nexp = mean_nexp,
      mean_pre_sd_nexp = mean_pre_sd_nexp, mean_sd_nexp = mean_sd_nexp,
      r_pre_post_exp = r_pre_post_exp, r_pre_post_nexp = r_pre_post_nexp,

      smd_to_cor = smd_to_cor, reverse_means_pre_post = reverse_means_pre_post,
      pre_post_to_smd = pre_post_to_smd
    )
  )
  es_means_se_pre_post <- with(
    x,
    es_from_means_se_pre_post(
      n_exp = n_exp, n_nexp = n_nexp,
      mean_pre_exp = mean_pre_exp, mean_exp = mean_exp,
      mean_pre_se_exp = mean_pre_se_exp, mean_se_exp = mean_se_exp,
      mean_pre_nexp = mean_pre_nexp, mean_nexp = mean_nexp,
      mean_pre_se_nexp = mean_pre_se_nexp, mean_se_nexp = mean_se_nexp,
      r_pre_post_exp = r_pre_post_exp, r_pre_post_nexp = r_pre_post_nexp,

      smd_to_cor = smd_to_cor, reverse_means_pre_post = reverse_means_pre_post,
      pre_post_to_smd = pre_post_to_smd
    )
  )
  es_means_ci_pre_post <- with(
    x,
    es_from_means_ci_pre_post(
      n_exp = n_exp, n_nexp = n_nexp,
      mean_pre_exp = mean_pre_exp, mean_exp = mean_exp,
      mean_pre_ci_lo_exp = mean_pre_ci_lo_exp, mean_pre_ci_up_exp = mean_pre_ci_up_exp,
      mean_ci_lo_exp = mean_ci_lo_exp, mean_ci_up_exp = mean_ci_up_exp,
      mean_pre_nexp = mean_pre_nexp, mean_nexp = mean_nexp,
      mean_pre_ci_lo_nexp = mean_pre_ci_lo_nexp, mean_pre_ci_up_nexp = mean_pre_ci_up_nexp,
      mean_ci_lo_nexp = mean_ci_lo_nexp, mean_ci_up_nexp = mean_ci_up_nexp,
      r_pre_post_exp = r_pre_post_exp, r_pre_post_nexp = r_pre_post_nexp,

      smd_to_cor = smd_to_cor, reverse_means_pre_post = reverse_means_pre_post,
      pre_post_to_smd = pre_post_to_smd
    )
  )

  es_means_change_sd <- with(x, es_from_mean_change_sd(
    n_exp = n_exp, n_nexp = n_nexp,
    mean_change_exp = mean_change_exp, mean_change_sd_exp = mean_change_sd_exp,
    mean_change_nexp = mean_change_nexp, mean_change_sd_nexp = mean_change_sd_nexp,
    r_pre_post_exp = r_pre_post_exp, r_pre_post_nexp = r_pre_post_nexp,

    smd_to_cor = smd_to_cor, reverse_mean_change = reverse_mean_change
  ))
  es_mean_change_se <- with(x, es_from_mean_change_se(
    n_exp = n_exp, n_nexp = n_nexp,
    mean_change_exp = mean_change_exp, mean_change_se_exp = mean_change_se_exp,
    mean_change_nexp = mean_change_nexp, mean_change_se_nexp = mean_change_se_nexp,
    r_pre_post_exp = r_pre_post_exp, r_pre_post_nexp = r_pre_post_nexp,

    smd_to_cor = smd_to_cor, reverse_mean_change = reverse_mean_change
  ))
  es_mean_change_ci <- with(x, es_from_mean_change_ci(
    n_exp = n_exp, n_nexp = n_nexp,
    mean_change_exp = mean_change_exp,
    mean_change_ci_lo_exp = mean_change_ci_lo_exp, mean_change_ci_up_exp = mean_change_ci_up_exp,
    mean_change_nexp = mean_change_nexp,
    mean_change_ci_lo_nexp = mean_change_ci_lo_nexp, mean_change_ci_up_nexp = mean_change_ci_up_nexp,
    r_pre_post_exp = r_pre_post_exp, r_pre_post_nexp = r_pre_post_nexp,

    smd_to_cor = smd_to_cor, reverse_mean_change = reverse_mean_change
  ))
  es_mean_change_pval <- with(x, es_from_mean_change_pval(
    n_exp = n_exp, n_nexp = n_nexp,
    mean_change_exp = mean_change_exp, mean_change_pval_exp = mean_change_pval_exp,
    mean_change_nexp = mean_change_nexp, mean_change_pval_nexp = mean_change_pval_nexp,
    r_pre_post_exp = r_pre_post_exp, r_pre_post_nexp = r_pre_post_nexp,

    smd_to_cor = smd_to_cor, reverse_mean_change = reverse_mean_change
  ))

  es_paired_t <- with(x, es_from_paired_t(
    paired_t_exp = paired_t_exp, paired_t_nexp = paired_t_nexp,
    n_exp = n_exp, n_nexp = n_nexp,
    r_pre_post_exp = r_pre_post_exp, r_pre_post_nexp = r_pre_post_nexp,

    smd_to_cor = smd_to_cor, reverse_paired_t = reverse_paired_t
  ))

  es_paired_t_pval <- with(x, es_from_paired_t_pval(
    paired_t_pval_exp = paired_t_pval_exp,
    paired_t_pval_nexp = paired_t_pval_nexp,
    n_exp = n_exp, n_nexp = n_nexp,
    r_pre_post_exp = r_pre_post_exp, r_pre_post_nexp = r_pre_post_nexp,

    smd_to_cor = smd_to_cor, reverse_paired_t_pval = reverse_paired_t_pval
  ))

  es_paired_f <- with(x, es_from_paired_f(paired_f_exp, paired_f_nexp,
    n_exp = n_exp, n_nexp = n_nexp,
    r_pre_post_exp = r_pre_post_exp, r_pre_post_nexp = r_pre_post_nexp,

    smd_to_cor = smd_to_cor, reverse_paired_f = reverse_paired_f
  ))

  es_paired_f_pval <- with(x, es_from_paired_f_pval(
    paired_f_pval_exp = paired_f_pval_exp,
    paired_f_pval_nexp = paired_f_pval_nexp,
    n_exp = n_exp, n_nexp = n_nexp,
    r_pre_post_exp = r_pre_post_exp, r_pre_post_nexp = r_pre_post_nexp,

    smd_to_cor = smd_to_cor, reverse_paired_f_pval = reverse_paired_f_pval
  ))
  # ANOVA, Student t-test  ----------------------------------------------
  es_t_student <- with(x, es_from_student_t(
    student_t = student_t, n_exp = n_exp, n_nexp = n_nexp,
    smd_to_cor = smd_to_cor, reverse_student_t = reverse_student_t
  ))
  es_t_student_pval <- with(x, es_from_student_t_pval(
    student_t_pval = student_t_pval, n_exp = n_exp, n_nexp = n_nexp,
    smd_to_cor = smd_to_cor, reverse_student_t_pval = reverse_student_t_pval
  ))
  es_anova_f <- with(x, es_from_anova_f(
    anova_f = anova_f, n_exp = n_exp, n_nexp = n_nexp,
    smd_to_cor = smd_to_cor, reverse_anova_f = reverse_anova_f
  ))
  es_anova_f_pval <- with(x, es_from_anova_pval(
    anova_f_pval = anova_f_pval, n_exp = n_exp, n_nexp = n_nexp,
    smd_to_cor = smd_to_cor, reverse_anova_f_pval = reverse_anova_f_pval
  ))
  es_etasq <- with(x, es_from_etasq(
    etasq = etasq, n_exp = n_exp, n_nexp = n_nexp,
    smd_to_cor = smd_to_cor, reverse_etasq = reverse_etasq
  ))
  es_etasq_adj <- with(x, es_from_etasq_adj(
    etasq_adj = etasq_adj, n_exp = n_exp, n_nexp = n_nexp, n_cov_ancova = n_cov_ancova,
    cov_outcome_r = cov_outcome_r, smd_to_cor = smd_to_cor, reverse_etasq = reverse_etasq
  ))

  # MD   ----------------------------------------------
  es_md_sd <- with(x, es_from_md_sd(
    md = md, md_sd = md_sd, n_exp = n_exp, n_nexp = n_nexp,
    smd_to_cor = smd_to_cor, reverse_md = reverse_md
  ))
  es_md_se <- with(x, es_from_md_se(
    md = md, md_se = md_se, n_exp = n_exp, n_nexp = n_nexp,
    smd_to_cor = smd_to_cor, reverse_md = reverse_md
  ))
  es_md_ci <- with(x, es_from_md_ci(
    md = md, md_ci_lo = md_ci_lo, md_ci_up = md_ci_up, n_exp = n_exp, n_nexp = n_nexp,
    smd_to_cor = smd_to_cor, reverse_md = reverse_md
  ))
  es_md_pval <- with(x, es_from_md_pval(
    md = md, md_pval = md_pval, n_exp = n_exp, n_nexp = n_nexp,
    smd_to_cor = smd_to_cor, reverse_md = reverse_md
  ))

  # RANGE/QUARTILES  ----------------------------------------------
  es_med_quarts <- with(x, es_from_med_quarts(
    q1_exp = q1_exp, med_exp = med_exp, q3_exp = q3_exp, n_exp = n_exp,
    q1_nexp = q1_nexp, med_nexp = med_nexp, q3_nexp = q3_nexp, n_nexp = n_nexp,
    smd_to_cor = smd_to_cor, reverse_med = reverse_med
  ))
  es_med_min_max <- with(x, es_from_med_min_max(
    min_exp = min_exp, med_exp = med_exp, max_exp = max_exp, n_exp = n_exp,
    min_nexp = min_nexp, med_nexp = med_nexp, max_nexp = max_nexp, n_nexp = n_nexp,
    smd_to_cor = smd_to_cor, reverse_med = reverse_med
  ))
  es_med_min_max_quarts <- with(x, es_from_med_min_max_quarts(
    min_exp = min_exp, q1_exp = q1_exp, med_exp = med_exp, q3_exp = q3_exp, max_exp = max_exp, n_exp = n_exp,
    min_nexp = min_nexp, q1_nexp = q1_nexp, med_nexp = med_nexp, q3_nexp = q3_nexp, max_nexp = max_nexp,
    n_nexp = n_nexp, smd_to_cor = smd_to_cor, reverse_med = reverse_med
  ))

  # ANCOVA MEANS ------------------------------------------------------------------
  es_ancova_means_sd <- with(x, es_from_ancova_means_sd(
    ancova_mean_exp = ancova_mean_exp, ancova_mean_nexp = ancova_mean_nexp,
    ancova_mean_sd_exp = ancova_mean_sd_exp, ancova_mean_sd_nexp = ancova_mean_sd_nexp,
    cov_outcome_r = cov_outcome_r, n_cov_ancova = n_cov_ancova, n_exp = n_exp, n_nexp = n_nexp,
    smd_to_cor = smd_to_cor, reverse_ancova_means = reverse_ancova_means
  ))
  es_ancova_means_se <- with(x, es_from_ancova_means_se(
    ancova_mean_exp = ancova_mean_exp, ancova_mean_nexp = ancova_mean_nexp,
    ancova_mean_se_exp = ancova_mean_se_exp, ancova_mean_se_nexp = ancova_mean_se_nexp,
    cov_outcome_r = cov_outcome_r, n_cov_ancova = n_cov_ancova, n_exp = n_exp, n_nexp = n_nexp,
    smd_to_cor = smd_to_cor, reverse_ancova_means = reverse_ancova_means
  ))
  es_ancova_means_ci <- with(x, es_from_ancova_means_ci(
    ancova_mean_exp = ancova_mean_exp, ancova_mean_nexp = ancova_mean_nexp,
    ancova_mean_ci_lo_exp = ancova_mean_ci_lo_exp, ancova_mean_ci_up_exp = ancova_mean_ci_up_exp,
    ancova_mean_ci_lo_nexp = ancova_mean_ci_lo_nexp, ancova_mean_ci_up_nexp = ancova_mean_ci_up_nexp,
    cov_outcome_r = cov_outcome_r, n_cov_ancova = n_cov_ancova, n_exp = n_exp, n_nexp = n_nexp,
    smd_to_cor = smd_to_cor, reverse_ancova_means = reverse_ancova_means
  ))
  es_ancova_means_sd_pooled_adj <- with(x, es_from_ancova_means_sd_pooled_adj(
    ancova_mean_exp = ancova_mean_exp, ancova_mean_nexp = ancova_mean_nexp,
    ancova_mean_sd_pooled = ancova_mean_sd_pooled,
    n_exp = n_exp, n_nexp = n_nexp, cov_outcome_r = cov_outcome_r, n_cov_ancova = n_cov_ancova,
    smd_to_cor = smd_to_cor, reverse_ancova_means = reverse_ancova_means
  ))
  es_ancova_means_sd_pooled <- with(x, es_from_ancova_means_sd_pooled_crude(
    ancova_mean_exp = ancova_mean_exp, ancova_mean_nexp = ancova_mean_nexp,
    mean_sd_pooled = mean_sd_pooled, cov_outcome_r = cov_outcome_r, n_cov_ancova = n_cov_ancova,
    n_exp = n_exp, n_nexp = n_nexp, smd_to_cor = smd_to_cor, reverse_ancova_means = reverse_ancova_means
  ))
  # plot
  es_plot_ancova_means <- with(x, es_from_plot_ancova_means(
    n_exp = n_exp, n_nexp = n_nexp,
    plot_ancova_mean_exp = plot_ancova_mean_exp, plot_ancova_mean_nexp = plot_ancova_mean_nexp,
    plot_ancova_mean_sd_lo_exp = plot_ancova_mean_sd_lo_exp, plot_ancova_mean_sd_lo_nexp = plot_ancova_mean_sd_lo_nexp,
    plot_ancova_mean_sd_up_exp = plot_ancova_mean_sd_up_exp, plot_ancova_mean_sd_up_nexp = plot_ancova_mean_sd_up_nexp,
    plot_ancova_mean_se_lo_exp = plot_ancova_mean_se_lo_exp, plot_ancova_mean_se_lo_nexp = plot_ancova_mean_se_lo_nexp,
    plot_ancova_mean_se_up_exp = plot_ancova_mean_se_up_exp, plot_ancova_mean_se_up_nexp = plot_ancova_mean_se_up_nexp,
    plot_ancova_mean_ci_lo_exp = plot_ancova_mean_ci_lo_exp, plot_ancova_mean_ci_lo_nexp = plot_ancova_mean_ci_lo_nexp,
    plot_ancova_mean_ci_up_exp = plot_ancova_mean_ci_up_exp, plot_ancova_mean_ci_up_nexp = plot_ancova_mean_ci_up_nexp,
    cov_outcome_r = cov_outcome_r, n_cov_ancova = n_cov_ancova,
    smd_to_cor = smd_to_cor, reverse_plot_ancova_means = reverse_plot_ancova_means
  ))
  # ANCOVA MD ------------------------------------------------------------------
  es_ancova_md_sd <- with(x, es_from_ancova_md_sd(ancova_md = ancova_md, ancova_md_sd = ancova_md_sd,
         cov_outcome_r = cov_outcome_r, n_cov_ancova = n_cov_ancova,
         n_exp = n_exp, n_nexp = n_nexp,
         smd_to_cor = smd_to_cor,
         reverse_ancova_md = reverse_ancova_md))
  es_ancova_md_se <- with(x, es_from_ancova_md_se(ancova_md = ancova_md, ancova_md_se = ancova_md_se,
         cov_outcome_r = cov_outcome_r, n_cov_ancova = n_cov_ancova,
         n_exp = n_exp, n_nexp = n_nexp,
         smd_to_cor = smd_to_cor,
         reverse_ancova_md = reverse_ancova_md))
  es_ancova_md_ci <- with(x, es_from_ancova_md_ci(ancova_md = ancova_md,
         ancova_md_ci_lo = ancova_md_ci_lo, ancova_md_ci_up = ancova_md_ci_up,
         cov_outcome_r = cov_outcome_r, n_cov_ancova = n_cov_ancova,
         n_exp = n_exp, n_nexp = n_nexp,
         smd_to_cor = smd_to_cor,
         reverse_ancova_md = reverse_ancova_md))
  es_ancova_md_pval <- with(x, es_from_ancova_md_pval(
        ancova_md = ancova_md,
        ancova_md_pval = ancova_md_pval,
        cov_outcome_r = cov_outcome_r, n_cov_ancova = n_cov_ancova,
        n_exp = n_exp, n_nexp = n_nexp,
        smd_to_cor = smd_to_cor,
        reverse_ancova_md = reverse_ancova_md))
  # ANCOVA F, T, P ------------------------------------------------------------------
  es_ancova_t <- with(x, es_from_ancova_t(
    ancova_t = ancova_t, cov_outcome_r = cov_outcome_r, n_cov_ancova = n_cov_ancova,
    n_exp = n_exp, n_nexp = n_nexp, smd_to_cor = smd_to_cor, reverse_ancova_t = reverse_ancova_t
  ))
  es_ancova_f <- with(x, es_from_ancova_f(
    ancova_f = ancova_f, cov_outcome_r = cov_outcome_r, n_cov_ancova = n_cov_ancova,
    n_exp = n_exp, n_nexp = n_nexp, smd_to_cor = smd_to_cor, reverse_ancova_f = reverse_ancova_f
  ))
  es_ancova_t_pval <- with(x, es_from_ancova_t_pval(
    ancova_t_pval = ancova_t_pval, cov_outcome_r = cov_outcome_r, n_cov_ancova = n_cov_ancova,
    n_exp = n_exp, n_nexp = n_nexp, smd_to_cor = smd_to_cor, reverse_ancova_t_pval = reverse_ancova_t_pval
  ))
  es_ancova_f_pval <- with(x, es_from_ancova_f_pval(
    ancova_f_pval = ancova_f_pval, cov_outcome_r = cov_outcome_r, n_cov_ancova = n_cov_ancova,
    n_exp = n_exp, n_nexp = n_nexp, smd_to_cor = smd_to_cor, reverse_ancova_f_pval = reverse_ancova_f_pval
  ))

  # CHI-SQ + PHI --------------------------------------------------------
  es_chisq <- with(x, es_from_chisq(
    chisq = chisq, n_sample = n_sample,
    n_cases = n_cases, n_exp = n_exp,
    yates_chisq = yates_chisq,
    reverse_chisq = reverse_chisq
  ))

  es_chisq_pval <- with(x, es_from_chisq_pval(
    chisq_pval = chisq_pval,
    n_cases = n_cases, n_exp = n_exp,
    yates_chisq = yates_chisq,
    n_sample = n_sample, reverse_chisq_pval = reverse_chisq_pval
  ))

  es_phi <- with(x, es_from_phi(
    phi = phi, n_sample = n_sample,
    n_cases = n_cases, n_exp = n_exp,
    reverse_phi = reverse_phi
  ))

  # COR-PB   --------------------------------------------------------
  es_r_point_bis <- with(x, es_from_pt_bis_r(
    pt_bis_r = pt_bis_r, n_exp = n_exp, n_nexp = n_nexp,
    smd_to_cor = smd_to_cor,
    reverse_pt_bis_r = reverse_pt_bis_r
  ))

  es_r_point_bis_pval <- with(x, es_from_pt_bis_r_pval(
    pt_bis_r_pval = pt_bis_r_pval,
    n_exp = n_exp, n_nexp = n_nexp,
    smd_to_cor = smd_to_cor,
    reverse_pt_bis_r_pval = reverse_pt_bis_r_pval
  ))

  # 2x2, PROP --------------------------------------------------------
  es_2x2 <- with(x, es_from_2x2(
    n_cases_exp = n_cases_exp, n_cases_nexp = n_cases_nexp,
    n_controls_exp = n_controls_exp, n_controls_nexp = n_controls_nexp,
    table_2x2_to_cor = table_2x2_to_cor, reverse_2x2 = reverse_2x2
  ))
  es_2x2_sum <- with(x, es_from_2x2_sum(
    n_cases_exp = n_cases_exp, n_cases_nexp = n_cases_nexp, n_exp = n_exp, n_nexp = n_nexp,
    table_2x2_to_cor = table_2x2_to_cor, reverse_2x2 = reverse_2x2
  ))
  es_prop <- with(x, es_from_2x2_prop(
    prop_cases_exp = prop_cases_exp, prop_cases_nexp = prop_cases_nexp, n_exp = n_exp, n_nexp = n_nexp,
    table_2x2_to_cor = table_2x2_to_cor, reverse_prop = reverse_prop
  ))

  # RR --------------------------------------------------------
  es_rr_se <- with(x, es_from_rr_se(
    rr = rr, logrr = logrr, logrr_se = logrr_se,
    baseline_risk = baseline_risk, n_exp = n_exp, n_nexp = n_nexp,
    n_cases = n_cases, n_controls = n_controls,
    rr_to_or = rr_to_or, smd_to_cor = smd_to_cor, reverse_rr = reverse_rr
  ))
  es_rr_ci <- with(x, es_from_rr_ci(
    rr = rr, logrr = logrr,
    baseline_risk = baseline_risk, n_exp = n_exp, n_nexp = n_nexp,
    rr_ci_lo = rr_ci_lo, logrr_ci_lo = logrr_ci_lo, rr_ci_up = rr_ci_up, logrr_ci_up = logrr_ci_up, n_cases = n_cases, n_controls = n_controls,
    rr_to_or = rr_to_or, smd_to_cor = smd_to_cor, reverse_rr = reverse_rr
  ))
  es_rr_pval <- with(x, es_from_rr_pval(
    rr = rr, logrr = logrr, rr_pval = rr_pval,
    baseline_risk = baseline_risk, n_exp = n_exp, n_nexp = n_nexp,
    n_cases = n_cases, n_controls = n_controls,
    rr_to_or = rr_to_or, smd_to_cor = smd_to_cor, reverse_rr = reverse_rr_pval
  ))

  # regression
  es_std_beta <- with(x, es_from_beta_std(
    beta_std = beta_std, sd_dv = sd_dv, n_exp = n_exp, n_nexp = n_nexp,
    smd_to_cor = smd_to_cor, reverse_beta_std = reverse_beta_std
  ))
  es_unstd_beta <- with(x, es_from_beta_unstd(
    beta_unstd = beta_unstd, sd_dv = sd_dv, n_exp = n_exp, n_nexp = n_nexp,
    smd_to_cor = smd_to_cor, reverse_beta_unstd = reverse_beta_unstd
  ))

  # user input
  es_user_crude <- es_from_user_crude(
    measure = measure,
    user_es_measure_crude = x$user_es_measure_crude,
    user_es_crude = x$user_es_crude,
    user_se_crude = x$user_se_crude,
    user_ci_lo_crude = x$user_ci_lo_crude,
    user_ci_up_crude = x$user_ci_up_crude
  )
  es_user_adj <- es_from_user_adj(
    user_es_measure_adj = x$user_es_measure_adj,
    user_es_adj = x$user_es_adj,
    user_se_adj = x$user_se_adj,
    user_ci_lo_adj = x$user_ci_lo_adj,
    user_ci_up_adj = x$user_ci_up_adj,
    measure = measure
  )
  # survival
  es_cases_time <- with(x, es_from_cases_time(
    n_cases_exp = n_cases_exp, n_cases_nexp = n_cases_nexp,
    time_exp = time_exp, time_nexp = time_nexp,
    reverse_irr = reverse_irr
  ))
  # variability
  var_means_sd <- with(x, es_variab_from_means_sd(
    mean_exp = mean_exp, mean_nexp = mean_nexp,
    mean_sd_exp = mean_sd_exp, mean_sd_nexp = mean_sd_nexp,
    n_exp = n_exp, n_nexp = n_nexp,
    reverse_means_variability = reverse_means_variability
  ))
  var_means_se <- with(x, es_variab_from_means_se(
    mean_exp = mean_exp, mean_nexp = mean_nexp,
    mean_se_exp = mean_se_exp, mean_se_nexp = mean_se_nexp,
    n_exp = n_exp, n_nexp = n_nexp,
    reverse_means_variability = reverse_means_variability
  ))
  var_means_ci <- with(x, es_variab_from_means_ci(
    mean_exp = mean_exp, mean_nexp = mean_nexp,
    mean_ci_lo_exp = mean_ci_lo_exp, mean_ci_up_exp = mean_ci_up_exp,
    mean_ci_lo_nexp = mean_ci_lo_nexp, mean_ci_up_nexp = mean_ci_up_nexp,
    n_exp = n_exp, n_nexp = n_nexp,
    reverse_means_variability = reverse_means_variability
  ))
  # ========================================= STEP 2. SYNTHESIS OF CALCULATIONS =========================================== #
  res <- list(
    es_user_crude = es_user_crude, es_user_adj = es_user_adj,
    # ES value
    es_cohen_d = es_cohen_d, es_cohen_d_adj = es_cohen_d_adj,
    es_hedges_g = es_hedges_g,
    es_odds_ratio = es_odds_ratio, es_odds_ratio_se = es_odds_ratio_se,
    es_odds_ratio_ci = es_odds_ratio_ci, es_odds_ratio_pval = es_odds_ratio_pval,
    es_pearson_r = es_pearson_r, es_fisher_z = es_fisher_z,
    # raw means + (SD SD_pooled SE CI)
    es_means_sd_raw = es_means_sd_raw, es_means_sd_pooled = es_means_sd_pooled,
    es_means_se_raw = es_means_se_raw, es_means_ci_raw = es_means_ci_raw,
    # paired MC
    es_means_change_sd = es_means_change_sd,
    es_mean_change_se = es_mean_change_se,
    es_mean_change_ci = es_mean_change_ci,
    es_mean_change_pval = es_mean_change_pval,
    # paired means
    es_means_sd_pre_post = es_means_sd_pre_post,
    es_means_se_pre_post = es_means_se_pre_post,
    es_means_ci_pre_post = es_means_ci_pre_post,
    # paired stats
    es_paired_t = es_paired_t, es_paired_t_pval = es_paired_t_pval,
    es_paired_f = es_paired_f, es_paired_f_pval = es_paired_f_pval,
    # t, z-score or F value
    es_t_student = es_t_student, es_t_student_pval = es_t_student_pval, es_anova_f = es_anova_f,
    es_anova_f_pval = es_anova_f_pval, es_etasq = es_etasq,
    # mean difference
    es_md_sd = es_md_sd, es_md_se = es_md_se, es_md_ci = es_md_ci, es_md_pval = es_md_pval,
    # raw medians + quarts/minmax
    es_med_quarts = es_med_quarts, es_med_min_max = es_med_min_max, es_med_min_max_quarts = es_med_min_max_quarts,
    # PLOT
    es_plot_means_raw = es_plot_means_raw,
    # ANCOVA means
    es_ancova_means_sd = es_ancova_means_sd,
    es_ancova_means_se = es_ancova_means_se,
    es_ancova_means_ci = es_ancova_means_ci,
    es_ancova_means_sd_pooled = es_ancova_means_sd_pooled,
    es_ancova_means_sd_pooled_adj = es_ancova_means_sd_pooled_adj,
    # ANCOVA stats
    es_ancova_t = es_ancova_t, es_ancova_f = es_ancova_f,
    es_ancova_t_pval = es_ancova_t_pval,
    es_ancova_f_pval = es_ancova_f_pval,
    es_etasq_adj = es_etasq_adj,
    es_plot_ancova_means = es_plot_ancova_means,
    # ANCOVA MD
    es_ancova_md_sd = es_ancova_md_sd,
    es_ancova_md_se = es_ancova_md_se,
    es_ancova_md_ci = es_ancova_md_ci,
    es_ancova_md_pval = es_ancova_md_pval,
    # CHI-SQ + PHI
    es_chisq = es_chisq,
    es_chisq_pval = es_chisq_pval,
    es_phi = es_phi, # es_phi_pval = es_phi_pval,
    # COR-PB
    es_r_point_bis = es_r_point_bis, es_r_point_bis_pval = es_r_point_bis_pval,
    # 2x2, PROP
    es_2x2 = es_2x2, es_2x2_sum = es_2x2_sum, es_prop = es_prop,
    # RR
    es_rr_se = es_rr_se, es_rr_ci = es_rr_ci, es_rr_pval = es_rr_pval,
    # regression
    es_std_beta = es_std_beta, es_unstd_beta = es_unstd_beta,
    # survival
    es_cases_time = es_cases_time,
    # variability
    var_means_sd = var_means_sd, var_means_se = var_means_se, var_means_ci = var_means_ci
    # var_means_change = var_means_change
  )

  class(res) <- "metaConvert"
  attr(res, "raw_data") <- x
  attr(res, "exp") <- exp
  attr(res, "measure") <- measure
  attr(res, "split_adjusted") <- split_adjusted
  attr(res, "es_selected") <- es_selected
  attr(res, "format_adjusted") <- format_adjusted
  attr(res, "hierarchy") <- hierarchy
  attr(res, "main_es") <- main_es
  return(res)
}
# x_save2 = x; list_df = df_es; ordering = ordering_crude; digits = digits;
# suffix = "_crude"; measure = measure
#
# x = dat; exp = TRUE; es_selected = "hierarchy"
# split_adjusted = TRUE; measure = "d_to_or, smd_to
# hierarchy = "means_sd"#hierarchy#"user_input_crude"
# digits = 3
# x$user_es_measure_crude = "ROR"
# x$user_es_measure_adj = "MOR"
# x$split_adjustedusted = x$se_adjusted = 2
# load_all(); View(convert_df(df.haza));


# x = dat
# measure = "d"
# main_es = TRUE
# split_adjusted = TRUE
# format_adjusted = "wide"
# verbose = TRUE
# es_selected = "hierarchy"
# phi_to_cor = "tetrachoric";
# chisq_to_cor = "tetrachoric";
# rr_to_or = "metaumbrella"
# or_to_rr = "metaumbrella_cases"
# or_to_cor = "bonett"
# cor_to_smd = "cooper"
# table_2x2_to_cor = "lipsey"
# yates_chisq = FALSE
# smd_to_cor = "viechtbauer"
# pre_post_to_smd = "morris"
# unit_type = "raw_scale"
# hierarchy = "rr_se"
# start_time <- Sys.time()
# end_time <- Sys.time()
# end_time - start_time

