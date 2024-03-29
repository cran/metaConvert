#' Aggregate a dataframe containing dependent effect sizes according to Borenstein's formulas
#'
#' @param x a dataframe that should be aggregated (must contain effect size values and standard errors).
#' @param dependence The type of dependence in your dataframe (can be either "outcomes" or "subgroups"). See details.
#' @param agg_fact A character string identifying the column name that contains the clustering units (all rows with the same \code{agg_fact} value will be aggregated together).
#' @param es A character string identifying the column name containing the effect size values. Default is "es".
#' @param se A character string identifying the column name containing the standard errors of the effect size. Default is "se".
#' @param cor_outcomes The correlation between effect sizes coming from the same clustering unit (only used when \code{dependence = "outcomes"}).
#' @param col_mean a vector of character strings identifying the column names for which the dependent values are summarized by taking their mean.
#' @param col_weighted_mean a vector of character strings identifying the column names for which the dependent values are summarized by taking their weighted mean.
#' @param weights The weights that will be used to estimated the weighted means.
#' @param col_sum a vector of character strings identifying the column names for which the dependent values are summarized by taking their sum.
#' @param col_min a vector of character strings identifying the column names for which the dependent values are summarized by taking their minimum.
#' @param col_max a vector of character strings identifying the column names for which the dependent values are summarized by taking their maximum.
#' @param col_fact a vector of character strings identifying the column names that are factors (different values will be separated by a "/" character).
#' @param na.rm a logical vector indicating whether missing values should be ignored in the calculations for the \code{col_mean}, \code{col_weighted_mean}, \code{col_sum}, \code{col_min} and \code{col_max} arguments.
#'
#' @details
#' 1. In the \code{dependence} argument, you should indicate "outcomes" if the dependence within the same clustering unit (e.g., study) is due to the presence of multiple effect sizes produced from the same participants (e.g., multiple outcomes, or multiple time-points)
#' 2. In the \code{dependence} argument, you should indicate "subgroups" if the dependence within the same clustering unit (e.g., study) is due to the presence of multiple effect sizes produced by independent subgroups (e.g., one effect size for boys, and one for girls).
#'
#' If you are working with ratio measures, make sure that the information on
#' the effect size estimates (i.e., the column passed to the es argument of the function)
#' is presented on the log scale.
#'
#' @return
#' The object returned by the \code{aggregate_df} contains, is a dataframe containing at the very least,
#' the aggregating factor, and the aggregated effect size values and standard errors. All columns indicated in the \code{col_*} arguments
#' will also be included in this dataframe.
#' \tabular{ll}{
#'  \code{row_id} \tab the row number in the original dataset.\cr
#'  \tab \cr
#'  \code{es} \tab the aggregated effect size value.\cr
#'  \tab \cr
#'  \code{se} \tab the standard error of the aggregated effect size.\cr
#'  \tab \cr
#'  \code{...} \tab any columns indicated in the \code{col_*} arguments.\cr
#'  \tab \cr
#' }
#'
#' @export aggregate_df
#'
#' @md
#'
#' @examples
#' res <- summary(convert_df(df.haza, measure = "d"))
#' aggregate_df(res, dependence = "outcomes", cor_outcomes = 0.8,
#'              agg_fact = "study_id", es = "es_crude", se = "se_crude",
#'              col_fact = c("outcome", "type_publication"))
#'
aggregate_df <- function(x, dependence = "outcomes", cor_outcomes = 0.8,
                         agg_fact, es = "es", se = "se",
                         col_mean = NA, col_weighted_mean = NA, weights = NA,
                         col_sum = NA, col_min = NA, col_max = NA,
                         col_fact = NA, na.rm = TRUE) {
  if (!dependence %in% c("outcomes", "subgroups")) {
    stop(paste0("'", dependence, "' not in tolerated measures. Possible inputs are: 'outcomes' or 'subgroups'."))
  }

  if ((missing(es) & !"es" %in% colnames(x)) | (missing(se) & !"se" %in% colnames(x))) {
    stop(paste0("Please indicate the column names storing information on the effect size values and standard errors using the 'es' & 'se' arguments"))
  }
  if ((missing(es) & "es" %in% colnames(x))) {
    es <- "es"
  }
  if ((missing(se) & "se" %in% colnames(x))) {
    se <- "se"
  }

  if (length(es) > 1 | length(se) > 1) {
    stop(paste0("The '", es, "' and the '", se, "' arguments should contain only 1 element"))
  }
  if (!es %in% colnames(x) | !se %in% colnames(x)) {
    stop(paste0("The '", es, "' or the '", se, "' variable has not been found in your dataset"))
  }

  if (dependence == "outcomes") {
    res <- .agg.outcomes(
      x = x, agg_fact = agg_fact, es = es, se = se,
      cor_outcomes = cor_outcomes,
      col_mean = col_mean,
      col_weighted_mean = col_weighted_mean, weights = weights,
      col_sum = col_sum, col_min = col_min, col_max = col_max,
      col_fact = col_fact, na.rm = na.rm
    )
  } else {
    res <- .agg.subgroups(
      x = x, agg_fact = agg_fact, es = es, se = se,
      col_mean = col_mean,
      col_weighted_mean = col_weighted_mean, weights = weights,
      col_sum = col_sum, col_min = col_min, col_max = col_max,
      col_fact = col_fact, na.rm = na.rm
    )
  }

  rownames(res) <- 1:nrow(res)
  return(res)
}

.agg.subgroups <- function(x, agg_fact, es, se,
                           col_mean = NA, col_weighted_mean = NA, weights = NA,
                           col_sum = NA, col_min = NA, col_max = NA,
                           col_fact = NA, na.rm = TRUE) {
  cols <- c(col_mean, col_weighted_mean, col_sum, col_min, col_fact)

  dup_col <- ifelse(length(na.omit(cols)) > 0, any(duplicated(na.omit(cols))), FALSE)

  if (dup_col) {
    stop("Each column of your dataset can be aggregated only once.")
  }

  if (any(!is.na(col_weighted_mean) & is.na(weights))) {
    stop("You should specify the 'weights' argument when requesting to calculate a weighted mean for some columns.")
  }

  df_mean <- df_w_mean <- df_sum <- df_min <- df_max <- df_fact <- NA

  dup_row <- duplicated(x[, agg_fact]) | duplicated(x[, agg_fact], fromLast = TRUE)

  x_tot <- cbind(x, es = x[, es], se = x[, se], row_index = 1:nrow(x))

  x_unique <- x_tot[!dup_row, ]

  x <- x_tot[dup_row, ]

  if (nrow(x) == length(unique(x[, agg_fact]))) {
    warning("the number of rows and unique 'agg_fact' values are equal")
    return(x_unique)
  }

  x_split <- split(x, x[, agg_fact])

  agg_bor <- do.call(rbind, lapply(x_split, .unique_es_subgroups)) # ,  measure = measure

  agg_bor$agg <- row.names(agg_bor)
  # agg_bor$agg <- unique(x[, agg_fact])

  if (any(!is.na(col_mean))) {
    df_mean <- aggregate.data.frame(x[, col_mean], by = list(agg = x[, agg_fact]), FUN = mean, na.rm = na.rm)
    colnames(df_mean) <- c("agg", col_mean)
  }
  if (any(!is.na(col_weighted_mean))) {
    x_pond <- x[, c(col_weighted_mean)] * x[, weights]

    sum_x_pond <- aggregate.data.frame(x_pond, by = list(agg = x[, agg_fact]), FUN = sum, na.rm = na.rm)
    sum_N <- aggregate.data.frame(x[, weights], by = list(agg = x[, agg_fact]), FUN = sum, na.rm = na.rm)

    df_w_mean <- cbind(unique(x[, agg_fact]), sum_x_pond$x / sum_N$x)
    colnames(df_w_mean) <- c("agg", col_weighted_mean)
  }

  if (any(!is.na(col_sum))) {
    df_sum <- aggregate.data.frame(x[, col_sum], by = list(agg = x[, agg_fact]), FUN = sum, na.rm = na.rm)
    colnames(df_sum) <- c("agg", col_sum)
  }

  if (any(!is.na(col_min))) {
    df_min <- aggregate.data.frame(x[, col_min], by = list(agg = x[, agg_fact]), FUN = min, na.rm = na.rm)
    colnames(df_min) <- c("agg", col_min)
  }

  if (any(!is.na(col_max))) {
    df_max <- aggregate.data.frame(x[, col_max], by = list(agg = x[, agg_fact]), FUN = max, na.rm = na.rm)
    colnames(df_max) <- c("agg", col_max)
  }

  if (any(!is.na(col_fact))) {
    df_fact <- aggregate.data.frame(x[, col_fact], by = list(agg = x[, agg_fact]), FUN = .agg.fact, na.rm = na.rm)
    colnames(df_fact) <- c("agg", col_fact)
  }

  df_add <- Reduce(
    function(x, y) merge(x, y),
    list(agg_bor, df_mean, df_w_mean, df_sum, df_min, df_max, df_fact)
  )

  df_add <- subset(df_add, select = -c(y))

  first <- c(
    which(colnames(df_add) == "row_index"),
    which(colnames(df_add) == "agg")
  )

  df_add_clean <- cbind(df_add[, first], df_add[, -first])

  colnames(df_add_clean)[colnames(df_add_clean) == "agg"] <- agg_fact

  df_return <- rbind(x_unique[, colnames(df_add_clean)], df_add_clean)

  res <- df_return[order(df_return$row_index), ]
  rownames(res) <- 1:nrow(res)

  return(res)
}

.agg.outcomes <- function(x, agg_fact, es, se,
                          col_mean, col_weighted_mean, weights,
                          col_sum, col_min, col_max,
                          col_fact, cor_outcomes, na.rm) {
  cols <- c(col_mean, col_weighted_mean, col_sum, col_min, col_fact)

  dup_col <- ifelse(length(na.omit(cols)) > 0, any(duplicated(na.omit(cols))), FALSE)

  if (dup_col) {
    stop("Each column of your dataset can be aggregated only once.")
  }

  if (any(!is.na(col_weighted_mean) & is.na(weights))) {
    stop("You should specify the 'weights' argument when requesting to calculate a weighted mean for some columns.")
  }

  df_mean <- df_w_mean <- df_sum <- df_min <- df_max <- df_fact <- NA


  dup_row <- duplicated(x[, agg_fact]) | duplicated(x[, agg_fact], fromLast = TRUE)

  x_tot <- cbind(x,
    es = x[, es], se = x[, se],
    # ci_lo = rep(NA_real_, nrow(x)), ci_up = rep(NA_real_, nrow(x)),
    row_index = 1:nrow(x)
  )

  x_unique <- x_tot[!dup_row, ]

  x <- x_tot[dup_row, ]

  if (nrow(x) == length(unique(x[, agg_fact]))) {
    warning("the number of rows and unique 'agg_fact' values are equal")
    return(x_unique)
  }

  x_split <- split(x, x[, agg_fact])
  # measure = "SMD"
  agg_bor <- do.call(rbind, lapply(x_split, .unique_es_outcomes, cor_outcomes = cor_outcomes)) # , measure = measure

  agg_bor$agg <- row.names(agg_bor)

  if (any(!is.na(col_mean))) {
    df_mean <- aggregate.data.frame(x[, col_mean], by = list(agg = x[, agg_fact]), FUN = mean, na.rm = na.rm)
    colnames(df_mean) <- c("agg", col_mean)
  }

  if (any(!is.na(col_weighted_mean))) {
    x_pond <- x[, c(col_weighted_mean)] * x[, weights]
    sum_x_pond <- aggregate.data.frame(x_pond, by = list(agg = x[, agg_fact]), FUN = sum, na.rm = na.rm)
    sum_N <- aggregate.data.frame(x[, weights], by = list(agg = x[, agg_fact]), FUN = sum, na.rm = na.rm)
    df_w_mean <- cbind(unique(x[, agg_fact]), sum_x_pond$x / sum_N$x)
    colnames(df_w_mean) <- c("agg", col_weighted_mean)
  }

  if (any(!is.na(col_sum))) {
    df_sum <- aggregate.data.frame(x[, col_sum], by = list(agg = x[, agg_fact]), FUN = sum, na.rm = na.rm)
    colnames(df_sum) <- c("agg", col_sum)
  }

  if (any(!is.na(col_min))) {
    df_min <- aggregate.data.frame(x[, col_min], by = list(agg = x[, agg_fact]), FUN = min, na.rm = na.rm)
    colnames(df_min) <- c("agg", col_min)
  }

  if (any(!is.na(col_max))) {
    df_max <- aggregate.data.frame(x[, col_max], by = list(agg = x[, agg_fact]), FUN = max, na.rm = na.rm)
    colnames(df_max) <- c("agg", col_max)
  }

  if (any(!is.na(col_fact))) {
    df_fact <- aggregate.data.frame(x[, col_fact], by = list(agg = x[, agg_fact]), FUN = .agg.fact, na.rm = na.rm)
    colnames(df_fact) <- c("agg", col_fact)
  }

  df_add <- Reduce(
    function(x, y) merge(x, y),
    list(
      agg_bor, df_mean,
      df_w_mean, df_sum,
      df_min, df_max, df_fact
    )
  )

  df_add <- subset(df_add, select = -c(y))

  first <- c(
    which(colnames(df_add) == "row_index"),
    which(colnames(df_add) == "agg")
  )

  df_add_clean <- cbind(df_add[, first], df_add[, -first])

  colnames(df_add_clean)[colnames(df_add_clean) == "agg"] <- agg_fact

  df_return <- rbind(x_unique[, colnames(df_add_clean)], df_add_clean)

  res <- df_return[order(df_return$row_index), ]
  return(res)
}

.unique_es_subgroups <- function(x) {
  weights <- 1 / (x$se^2)
  es_list <- x$es
  # if (measure == "SMD") { #, measure
  #   es_list = x$es
  # } else {
  #   es_list = log(x$es)
  # }

  mean_es <- sum(weights * es_list) / sum(weights)
  se <- sqrt(1 / sum(weights))
  # es = ifelse(TRUE, mean_es, exp(mean_es))
  es <- mean_es

  res <- data.frame(
    es = es,
    se = se,
    row_index = x$row_index[1]
  )

  return(res)
}
.unique_es_outcomes <- function(x, cor_outcomes) { # , measure
  var_es <- x$se^2
  prod_se <- x$se %*% t(x$se)
  prod_se_r <- prod_se * cor_outcomes
  prod_se_r[lower.tri(prod_se_r)] <- 0
  diag(prod_se_r) <- 0

  # mean_es <- ifelse(TRUE, mean(x$es), exp(mean(log(x$es))))
  mean_es <- mean(x$es)
  var <- (1 / length(x$es))^2 * (sum(var_es) + 2 * sum(prod_se_r))
  se <- sqrt(var)


  res <- data.frame(
    es = mean_es,
    se = se,
    row_index = x$row_index[1]
  )

  return(res)
}

.agg.fact <- function(x, na.rm = na.rm) {
  if (na.rm) {
    ifelse(length(unique(na.omit(x))) == 0,
      "N/A",
      ifelse(length(unique(na.omit(x))) == 1,
        unique(na.omit(x)),
        paste(na.omit(sort(unique(x))), collapse = " / ")
      )
    )
  } else {
    ifelse(length(unique(x)) == 0,
      "N/A",
      ifelse(length(unique(x)) == 1,
        unique(x),
        paste(sort(unique(x)), collapse = " / ")
      )
    )
  }
}
