#' Statistical helper functions for proteomics analysis workflows
#'
#' This file contains small utility functions used throughout the spheroid
#' proteomics analysis pipeline. The helpers support outlier detection,
#' differential expression testing, feature filtering, and principal component
#' analysis for `SummarizedExperiment` objects.
#'
#' The functions included in this file are designed to streamline exploratory and
#' statistical analysis of proteomics data and are intended to be used in
#' downstream plotting and reporting workflows.
#'
#' @name stats_funcs
#' @aliases stats_funcs
#' @title Statistical helper functions for proteomics datasets
#' @description Utility functions for outlier detection, differential expression
#' analysis, feature-level classification, and PCA on SummarizedExperiment
#' objects.
#' @keywords internal

library(dplyr)
library(magrittr)
library(tibble)
library(tidyr)
library(SummarizedExperiment)


#' Identify z-score outliers in a numeric vector
#'
#' Computes a z-score for each value in `x` and flags observations that exceed a
#' specified sensitivity threshold. This is useful for quickly identifying values
#' that deviate substantially from the center of a distribution.
#'
#' @param x Numeric vector whose values should be evaluated for outliers.
#' @param sensitivity Numeric threshold used to define how extreme a z-score
#'   must be to be considered an outlier. Defaults to 3.
#' @param sided Character string indicating whether to flag two-sided,
#'   high-end-only, or low-end-only outliers. Must be one of `"two"`, `"right"`,
#'   or `"left"`.
#'
#' @return A logical vector the same length as `x`, where `TRUE` indicates an
#'   outlier according to the chosen detection direction and sensitivity.
#'
#' @details The function standardizes each value by subtracting the mean and
#' dividing by the standard deviation, both computed with `na.rm = TRUE`. A
#' two-sided test flags values with $|z| > sensitivity$, while directional tests
#' flag only high or low tail values.
#'
#' @examples
#' values <- c(10, 11, 12, 13, 100)
#' get_zscore_outliers(values)
#' get_zscore_outliers(values, sensitivity = 2, sided = "right")
#' get_zscore_outliers(values, sided = "left")
get_zscore_outliers <- function(x, sensitivity = 3, sided = "two") {
  z_scores <- (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  if (sided == "two") {
    return(abs(z_scores) > sensitivity)
  } else if (sided == "right") {
    return(z_scores > sensitivity)
  } else if (sided == "left") {
    return(z_scores < -sensitivity)
  } else {
    stop("Invalid value for 'sided'. Choose from 'two', 'right', 'left'.")
  }
}


#' Differential expression testing between two sample groups
#'
#' Runs a limma-based comparison between two groups in a
#' \code{SummarizedExperiment} and appends the estimated statistics back to the
#' feature metadata in \code{rowData(experiment)}. The comparison is performed
#' on the assay specified by \code{assay_name}, with samples assigned to the
#' design matrix in the order required by limma.
#'
#' @param experiment A \code{SummarizedExperiment} object containing expression
#'   values in one or more assays and sample group labels in
#'   \code{colData(experiment)$group}.
#' @param group_1 Character string identifying the reference group in
#'   \code{colData(experiment)[["group"]]}.
#' @param group_2 Character string identifying the comparison group in
#'   \code{colData(experiment)[["group"]]}.
#' @param assay_name Character string naming the assay to analyze. Defaults to
#'   \code{"proteomics"}.
#' @param align_by Column name in \code{rowData(experiment)} used to align the
#'   limma output back to the feature metadata. Defaults to
#'   \code{"Protein.Group"}.
#'
#' @return The original \code{experiment} object with limma statistics merged
#'   into \code{rowData(experiment)} using \code{align_by}. The result table
#'   includes all columns returned by \code{limma::topTable()} plus the feature
#'   identifier column used for alignment.
#'
#' @details
#' The function identifies sample indices for the two specified groups, builds a
#' two-group design matrix, fits the model with limma, and adds the resulting
#' differential expression table back to the feature metadata. A stop condition
#' is raised if \code{align_by} is not found among the columns of
#' \code{rowData(experiment)}.
#'
#' @examples
#' \dontrun{
#' # se is a SummarizedExperiment with each sample labeled by group in
#' # colData(se)$group and a proteomics assay.
#' se <- calc_diff_ex(
#'   experiment = se,
#'   group_1 = "control",
#'   group_2 = "treatment",
#'   assay_name = "proteomics",
#'   align_by = "Protein.Group"
#' )
#'
#' # Inspect the added statistics
#' rowData(se) |> as_tibble() |> select(Protein.Group, logFC, adj.P.Val)
#' }
calc_diff_ex <- function(
  experiment,
  group_1,
  group_2,
  assay_name = "proteomics",
  align_by = "Protein.Group"
) {
  # making sure that the 'align_by' group is present in the SummarizedExperiment
  if (!align_by %in% colnames(rowData(experiment))) {
    stop(
      base::paste(
        align_by,
        "not in rowData. Available columns:",
        paste0(colnames(rowData(experiment)), collapse = ","),
        sep = " "
      ),
      call. = FALSE
    )
  }

  # retrieving the samples (column ids) of the groups we want to use for the
  # differential expression analysis
  samples_1 <- which(
    SummarizedExperiment::colData(experiment)[["group"]] == group_1
  )
  samples_2 <- which(
    SummarizedExperiment::colData(experiment)[["group"]] == group_2
  )

  # next we set up the design matrix for the fit function used by limma to
  # determine the differentially expressed features. The design matrix uses
  # a factor based on the samples
  fac <- factor(rep(c(2, 1), c(length(samples_2), length(samples_1))))
  design <- stats::model.matrix(~fac)

  # next we extract the assay of interest from the `SummarizedExperiment`
  # object and create the differential expression result table. Note that we
  # switch the order of groups here for `limma`
  dat <- SummarizedExperiment::assays(experiment)[[assay_name]]

  result_table <- dat[, c(samples_2, samples_1)] %>%
    limma::lmFit(., design) %>%
    limma::eBayes(.) %>%
    limma::topTable(., coef = 2, number = Inf, sort.by = "none")

  # The assumption is that the input `SummarizedExperiment` had only one assay
  # associated with it, hence we will just augment the object and return it
  # containing the results from limma.

  result_table <- result_table %>%
    rownames_to_column(., var = align_by) %>%
    as_tibble(.)

  # This makes sure that the differential expression results for each feature
  # are properly aligned

  SummarizedExperiment::rowData(experiment) <- left_join(
    as_tibble(SummarizedExperiment::rowData(experiment)),
    result_table,
    by = align_by
  )

  return(experiment)
}


#' Extract and classify differentially expressed features
#'
#' Summarizes differential expression results stored in `rowData(diff_ex)` by
#' selecting a feature identifier column plus the log fold-change and adjusted
#' p-value columns, then classifying each feature as up-regulated, down-
#' regulated, or not significant based on thresholds.
#'
#' @param diff_ex A `SummarizedExperiment` object whose `rowData()` contains the
#'   differential expression statistics.
#' @param feature_type Character string naming the column in `rowData(diff_ex)`
#'   that uniquely identifies each feature.
#' @param log_fc Numeric threshold for the absolute log fold-change used to
#'   declare a feature differentially expressed. Defaults to `1`.
#' @param pval Numeric adjusted p-value threshold for significance. Defaults to
#'   `0.05`.
#'
#' @return A tibble with one row per feature and columns for the feature label,
#'   `logFC`, `adj.P.Val`, `expression_change`, and `labels`. Features meeting
#'   the fold-change and significance criteria are labeled as `"up"` or
#'   `"down"`; all others are labeled as `"not significant"`.
#'
#' @details The function validates that `feature_type` exists in the feature
#' metadata, then selects the required columns and adds a derived
#' `expression_change` field using `case_when()`. The `labels` column is filled
#' only for features with a significant change, making it useful for annotation
#' or plotting.
#'
#' @examples
#' \dontrun{
#' # diff_ex is a SummarizedExperiment with rowData columns including
#' # "Protein.Group", "logFC", and "adj.P.Val".
#' results <- extract_diff_ex_features(
#'   diff_ex = diff_ex,
#'   feature_type = "Protein.Group",
#'   log_fc = 1,
#'   pval = 0.05
#' )
#'
#' head(results)
#' }
extract_diff_ex_features <- function(
  diff_ex,
  feature_type,
  log_fc = 1,
  pval = 0.05
) {
  df_ret <- rowData(diff_ex) |> tibble::as_tibble()

  if (!feature_type %in% names(df_ret)) {
    stop(paste0(
      "feature_type '",
      feature_type,
      "' not present in rowData of diff_ex",
    ))
  }

  df_ret %<>%
    dplyr::select(
      feature = feature_type,
      'logFC',
      'adj.P.Val'
    ) %>%
    mutate(
      expression_change = case_when(
        (logFC >= log_fc & adj.P.Val < pval) ~ 'up',
        (logFC <= -log_fc & adj.P.Val < pval) ~ 'down',
        .default = 'not significant',
      )
    ) %>%
    mutate(
      labels = case_when(
        expression_change %in% c('up', 'down') ~ feature,
        .default = NA
      )
    )

  return(df_ret)
}


#' Generate principal components for a SummarizedExperiment
#'
#' Runs a principal component analysis on the primary assay of a
#' `SummarizedExperiment` and returns the fitted model along with sample-level
#' component scores and variance explained by each principal component.
#'
#' @param summarized_experiment A `SummarizedExperiment` object containing
#'   expression values in one or more assays.
#'
#' @return A named list with the following elements:
#' \itemize{
#'   \item{`prcomp_res`}{The full result object returned by `stats::prcomp()`.}
#'   \item{`prcs`}{A data frame of principal component scores for each sample,
#'     combined with the sample metadata in `colData(summarized_experiment)`.}
#'   \item{`explained_variance`}{A one-row data frame containing the proportion
#'     of variance explained by each principal component.}
#' }
#'
#' @details The function extracts the first assay from the input object,
#' transposes it so that samples are rows and features are columns, then calls
#' `prcomp()`. The resulting component scores are combined with sample metadata
#' for downstream plotting or clustering analyses.
#'
#' @examples
#' \dontrun{
#' pca_results <- generate_pca(se)
#' str(pca_results)
#' head(pca_results$prcs)
#' }
generate_pca <- function(summarized_experiment) {
  prc_res <- summarized_experiment |>
    assay() |> # extracting the `assay` (e.g. proteomics) from the SE
    t() |> # transposing `assay` object to fit with prcomp requirements
    prcomp() %$% # prcomp requires samples = rows; features = columns
    list(
      prcomp_res = .,
      prcs = cbind(x, colData(summarized_experiment)),
      explained_variance = as.data.frame(summary(.)$importance)[2, ]
    )

  return(prc_res)
}
