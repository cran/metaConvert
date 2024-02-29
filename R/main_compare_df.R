#' Flag the differences between two dataframes.
#'
#' @param df_extractor_1 a first dataset. Differences with the second dataset will be flagged in green.
#' @param df_extractor_2 a second dataset. Differences with the first dataset will be flagged in red.
#' @param tolerance the cut-off value used to flag differences between two numeric values
#' @param tolerance_type must be either 'difference' or 'ratio'
#' @param output type of object returned by the function (see 'Value' section). Must be either 'wide', 'long', 'html', 'html2' or 'xlsx'.
#' @param file_name the name of the generated file (only used when \code{output="xlsx"})
#' @param ordering_columns column names that should be used to re-order the two datasets before running the comparisons
#'
#' @details
#' This function aims to facilitate the comparison of two datasets created by blind data extractors during a systematic review.
#' It is a wrapper of several functions from the 'compareDF' package.
#'
#' @return
#' This function returns a dataframe composed of the rows that include a
#' difference (all identical rows are removed).
#' Several outputs can be requested :
#' 1. setting \code{output="xlsx"} returns an excel file. A message indicates the location of the generated file on your computer.
#' 1. setting \code{output="html"} returns an html file
#' 1. setting \code{output="html2"} returns an html file (only useful when the "html" command did not make the html pane appear in R studio).
#' 1. setting \code{output="wide"} a wide dataframe
#' 1. setting \code{output="long"} a long dataframe
#'
#' @md
#'
#' @export compare_df
#'
#' @references
#' Alex Joseph (2022). compareDF: Do a Git Style Diff of the Rows Between Two Dataframes with Similar Structure. R package version 2.3.3. https://CRAN.R-project.org/package=compareDF
#'
#' @examples
#' df.compare1 = df.compare1[order(df.compare1$author), ]
#' df.compare2 = df.compare2[order(df.compare2$year), ]
#' names(df.compare1)[2] <- "generate_warning"
#'
#' compare_df(
#'   df_extractor_1 = df.compare1,
#'   df_extractor_2 = df.compare2,
#'   ordering_columns = c("study_id")
#' )

compare_df <- function(df_extractor_1, df_extractor_2,
                       ordering_columns = NULL,
                       tolerance = 0,
                       tolerance_type = "ratio",
                       output = "html",
                       file_name = "comparison.xlsx") {

  df_extractor_1 = data.frame(df_extractor_1)
  df_extractor_2 = data.frame(df_extractor_2)

  cols_available = unique(append(names(df_extractor_1), names(df_extractor_2)))
  cols_retained = intersect(names(df_extractor_1), names(df_extractor_2))
  cols_unique = cols_available[!cols_available %in% cols_retained]
  if (length(cols_unique) > 0) {
    warning("\nThe columns ", paste(cols_unique, collapse = " / "),
            " were removed from the comparison since they were included only in one of the two datasets.\n")
  }


  if (!is.null(ordering_columns)) {
    df_extractor_1 = suppressWarnings(df_extractor_1[order(df_extractor_1[, ordering_columns]), cols_retained])
    df_extractor_2 = suppressWarnings(df_extractor_2[order(df_extractor_2[, ordering_columns]), cols_retained])
    } else {
    df_extractor_1 = df_extractor_1[, cols_retained]
    df_extractor_2 = df_extractor_2[, cols_retained]
  }
  df_extractor_1$rowname <- 1:nrow(df_extractor_1)
  df_extractor_2$rowname <- 1:nrow(df_extractor_2)

  comp <- suppressMessages(compareDF::compare_df(
    df_new = df_extractor_1,
    df_old = df_extractor_2,
    group_col = "rowname",
    exclude = NULL,
    tolerance = tolerance,
    tolerance_type = tolerance_type,
    stop_on_error = TRUE,
    change_markers = c("df_extractor_1", "df_extractor_2", "="),
    round_output_to = 3
  ))

  if (output == "long") {
    res <- comp$comparison_df[order(as.numeric(as.character(comp$comparison_df$rowname))), ]
  } else if (output == "wide") {
    res <- compareDF::create_wide_output(comp)
    res <- res[order(as.numeric(as.character(res$rowname))), ]
  } else if (output == "html") {
    res <- suppressMessages(compareDF::create_output_table(comp, output_type = "html", limit = 5000))
  } else if (output == "html2") {
    res <- suppressMessages(compareDF::view_html(comp))
  } else if (output == "xlsx") {
    if (!grepl(".xlsx", file_name)) {
      file_name <- paste0(file_name, ".xlsx")
    }
    message("File generated at '", paste0(getwd(), "/", file_name, "'"))
    res <- compareDF::create_output_table(comp,
      output_type = "xlsx",
      file_name = file_name, limit = 5000
    )
  }
  return(res)
}
