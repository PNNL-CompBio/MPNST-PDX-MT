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


extract_diff_ex_features <- function(
  diff_ex,
  feature_type,
  log_fc = 1,
  pval = 0.05
) {
  df_ret <- rowData(diff_ex) |> as_tibble()

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

  df_ret
}


#' Generate principal components for SummarizedExperiment
#'
#' Helper function to generate principal components for an SummarizedExperiment
#' object.
#'
#' @param summarized_experiment The {SummarizedExperiment} object for which the
#'   principle components should be generated
#' @returns
#' Returns a named list containing the following components
#'
#' * prcomp_res The full output of `stats.prcomp()`
#' * prcs A `data.frame` containing the combined principal component values per
#'   sample including additional sample data extracted from
#'   `colData(summarized_experiment)`
#' * explained_variance A `data.frame` containing the portion of variance for
#'   each PC
#'
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

  prc_res
}
