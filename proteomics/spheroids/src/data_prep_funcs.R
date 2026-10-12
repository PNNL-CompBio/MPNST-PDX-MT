library(tibble)
library(dplyr)
library(tidyr)
library(magrittr)

library(SummarizedExperiment)
library(limma)
library(matrixStats)
library(sva)


prepare_proteomics_data <- function(
  proteomics_data,
  proteomics_meta,
  type = c('global', 'phospho'),
  ...
) {
  # checking that either global or phospho is selected for type
  # match.arg() will throw an error otherwise and execution will stop
  sel_type <- match.arg(type)

  # Handeling extra arguments passed to the function.
  # Optional additional arguments are samp_cutoff, missingness_cutoff, and group_contains_spheroid
  extra_args <- list(...)

  # Checking for sample_col_prefix argument in extra_args.
  # Needed to extract the sample columns from the proteomics data.
  # If it is not present, we will use the default value of "cNF".
  # TODO: Make mandatory?
  if (is.null(extra_args$sample_col_prefix)) {
    sample_col_prefix <- "cNF"
  } else {
    sample_col_prefix <- extra_args$sample_col_prefix
  }

  if (
    !is.null(extra_args$sample_missingness_cutoff) &&
      !is.null(extra_args$samples_to_exclude)
  ) {
    stop(
      "Cannot specify both `sample_missingness_cutoff` and `samples_to_exclude`.",
      call. = FALSE
    )
  }

  if (!is.null(extra_args$sample_missingness_cutoff)) {
    sample_missingness_cutoff <- extra_args$sample_missingness_cutoff
    if (sample_missingness_cutoff < 0 | sample_missingness_cutoff > 1) {
      stop(
        base::paste(
          "sample_missingness_cutoff needs to be '0 <= sample_missingness_cutoff <=1'",
          "or 'NULL'. Supplied value:",
          sample_missingness_cutoff
        ),
        call. = FALSE
      )
    }
  } else {
    sample_missingness_cutoff <- NULL
  }

  if (!is.null(extra_args$samples_to_exclude)) {
    samples_to_exclude <- extra_args$samples_to_exclude
    if (!is.character(samples_to_exclude)) {
      stop(
        base::paste(
          "samples_to_exclude needs to be a character vector of sample names",
          "or 'NULL'. Supplied value:",
          samples_to_exclude
        ),
        call. = FALSE
      )
    }
  } else {
    samples_to_exclude <- NULL
  }

  if (!is.null(extra_args$feature_missingness_cutoff)) {
    feature_missingness_cutoff <- extra_args$feature_missingness_cutoff
    if (feature_missingness_cutoff < 0 | feature_missingness_cutoff > 1) {
      stop(
        base::paste(
          "feature_missingness_cutoff needs to be '0 <= feature_missingness_cutoff <=1'",
          "or 'NULL'. Supplied value:",
          feature_missingness_cutoff
        ),
        call. = FALSE
      )
    }
  } else {
    feature_missingness_cutoff <- NULL
  }

  if (!is.null(extra_args$imputation_method)) {
    imputation_method <- match.arg(
      extra_args$imputation_method,
      c('SampMin', 'means', 'median')
    )
  } else {
    imputation_method <- NULL
  }

  if (!is.null(extra_args$group_contains_spheroid)) {
    group_contains_spheroid <- extra_args$group_contains_spheroid
    if (
      !is.logical(group_contains_spheroid) ||
        length(group_contains_spheroid) != 1L
    ) {
      stop("`group_contains_spheroid` must be a single TRUE or FALSE value.")
    }
  }

  # cleaning up the metadata by splitting the condition column into group, replicate, and if applicable spheroid, timepoint, and agent
  proteomics_meta_clean <- clean_meta_data(
    data = proteomics_meta,
    group_contains_spheroid = group_contains_spheroid
  )

  if (sel_type == 'global') {
    proteomics_data_cleaned <- clean_global_proteomics_data(
      data = proteomics_data,
      sample_col_prefix = sample_col_prefix
    )
  } else if (sel_type == 'phospho') {
    proteomics_data_cleaned <- clean_phospho_proteomics_data(
      data = proteomics_data,
      sample_col_prefix = sample_col_prefix
    )
  }
  proteomics_data_clean <- proteomics_data_cleaned$data
  proteomicsMapping <- proteomics_data_cleaned$mapping

  experiment <- SummarizedExperiment(
    assays = list(proteomics = proteomics_data_clean),
    rowData = proteomicsMapping,
    colData = proteomics_meta_clean[colnames(proteomics_data_clean), ]
  )

  if (!is.null(feature_missingness_cutoff)) {
    experiment %<>%
      filter_for_feature_missingness(
        .,
        feature_missingness_cutoff
      )
  }

  if (!is.null(samples_to_exclude)) {
    experiment <- experiment[, !colnames(experiment) %in% samples_to_exclude]
  }

  if (!is.null(sample_missingness_cutoff)) {
    experiment %<>%
      filter_for_sample_cutoff(
        .,
        sample_missingness_cutoff
      )
  }

  if (is.null(imputation_method)) {
    experiment %<>%
      log_transform(., assay = "proteomics")
  } else {
    experiment %<>%
      impute_proteomics(., assay = "proteomics", imputation_method) %>%
      log_transform(., assay = "proteomics")
  }

  return(experiment)
}


clean_meta_data <- function(
  data,
  group_contains_spheroid = FALSE
) {
  if (group_contains_spheroid == TRUE) {
    data_cleaned <- data %>%
      separate_wider_regex(
        .,
        col = "condition",
        c("group" = ".*", "_", "replicate" = ".*"),
        cols_remove = FALSE
      ) %>%
      separate_wider_delim(
        .,
        col = "group",
        delim = "_",
        names = c("spheroid", "timepoint", "agent"),
        too_few = "align_start",
        cols_remove = FALSE
      ) |>
      tibble::column_to_rownames('sample')
  } else {
    data_cleaned <- data %>%
      separate_wider_regex(
        .,
        col = "condition",
        c("group" = ".*", "_", "replicate" = ".*"),
        cols_remove = FALSE
      ) %>%
      separate_wider_delim(
        .,
        col = "group",
        delim = "_",
        names = c("timepoint", "agent"),
        too_few = "align_start",
        cols_remove = FALSE
      ) |>
      tibble::column_to_rownames('sample')
  }
  return(data_cleaned)
}

clean_global_proteomics_data <- function(
  data,
  sample_col_prefix = "cNF"
) {
  prot_data_cleaned <- data %>%
    column_to_rownames(., var = "Protein.Group") %>%
    select(., contains(sample_col_prefix)) |>
    rownames_to_column("Protein.Group")

  prot_mapping <- data[0:4] %>%
    separate_wider_delim(
      .,
      col = "Genes",
      names = c("Genes", NA),
      delim = ";",
      cols_remove = TRUE,
      too_few = "align_start",
      too_many = "drop"
    ) %>%
    semi_join(., prot_data_cleaned, by = "Protein.Group")

  prot_data_cleaned %<>% column_to_rownames(., var = "Protein.Group")
  return(list(
    'data' = prot_data_cleaned,
    'mapping' = prot_mapping
  ))
}

clean_phospho_proteomics_data <- function(
  data,
  sample_col_prefix = "cNF"
) {
  prot_data_cleaned <- data %>%
    unite(., "Phospho.Site", Residue, Site, sep = "", remove = TRUE) %>%
    unite(
      .,
      "Protein.With.Phospho.Site",
      Protein,
      Phospho.Site,
      sep = "-"
    ) %>%
    column_to_rownames(., var = "Protein.With.Phospho.Site") %>%
    within(., rm("Protein.Names", "Gene.Names", "Sequence")) %>%
    rownames_to_column('Protein.With.Phospho.Site')

  prot_mapping <- data[0:5] %>%
    separate_wider_delim(
      .,
      col = "Gene.Names",
      names = c("Gene.Names", NA),
      delim = ";",
      cols_remove = TRUE,
      too_few = "align_start",
      too_many = "drop"
    ) %>%
    unite(., "Phospho.Site", Residue, Site, sep = "", remove = TRUE) %>%
    unite(
      .,
      "Protein.With.Phospho.Site",
      Protein,
      Phospho.Site,
      sep = "-",
      remove = FALSE
    ) %>%
    unite(
      .,
      "Gene.With.Phospho.Site",
      Gene.Names,
      Phospho.Site,
      sep = "-",
      remove = FALSE
    ) %>%
    semi_join(., prot_data_cleaned, by = "Protein.With.Phospho.Site")

  prot_data_cleaned %<>%
    column_to_rownames(., var = "Protein.With.Phospho.Site")

  #calculate median abundnace for each phosphosite and rank
  medAbund <- apply(prot_data_cleaned, 1, function(x) {
    median(x, na.rm = T)
  })
  qvals <- quantile(medAbund, probs = c(0.1, 0.2, 0.5, 1), na.rm = TRUE)
  siteQuantile <- sapply(medAbund, function(x) {
    names(qvals)[which(x <= qvals)[1]]
  })

  prot_mapping$medAbundance <- medAbund
  prot_mapping$quantile <- siteQuantile
  return(list(
    'data' = prot_data_cleaned,
    'mapping' = prot_mapping
  ))
}


#' Log-transform a SummarizedExperiment assay
#'
#' Applies a base-2 log transformation to the selected assay in a
#' `SummarizedExperiment` object and returns the updated object.
#'
#' @param experiment A `SummarizedExperiment` object that contains the assay to
#'   transform.
#' @param assay A character scalar naming the assay to log-transform.
#'   Defaults to `"proteomics"`.
#' @returns A `SummarizedExperiment` object with the specified assay replaced by
#'   its log2-transformed values.
#'
log_transform <- function(experiment, assay = "proteomics") {
  proteomics_data <- SummarizedExperiment::assay(experiment, assay) %>%
    mutate(across(everything(), as.numeric)) %>% # turn everything numeric
    log2(.) # log transform all cells

  SummarizedExperiment::assay(
    experiment,
    withDimnames = FALSE
  ) <- proteomics_data
  SummarizedExperiment::assayNames(experiment) <- assay
  return(experiment)
}

impute_proteomics <- function(
  experiment,
  assay = 'proteomics',
  imputation_method = 'SampMin'
) {
  # retrieve assay from SummarizedExperiment and make sure that all values
  # are numeric such that we can impute
  proteomics_data <- SummarizedExperiment::assays(experiment)[[assay]] %>%
    mutate(across(everything(), as.numeric)) # turn every cell numeric

  if (imputation_method == "SampMin") {
    # first calculates the minimum protein intensity per sample (column) using
    # `matrixStats::colMins()`
    sample_mins <- matrixStats::colMins(
      as.matrix(proteomics_data),
      na.rm = TRUE
    )
    # then NAs in each sample (column) to the sample minimum (stored as list)
    proteomics_data <- proteomics_data %>%
      replace_na(., as.list(sample_mins))
  } else if (imputation_method == "means") {
    # calculates the mean abundance of a given feature (protein) across all
    # samples (columns) via `base::rowMeans()` stores it temporarily into a new
    # column `imputeVal` and then iterates over all columns (excluding the newly
    # generated column) to set N/As to the imputation value using
    # `dplyr::coalesce`. Finally the temporary column is removed from the data.
    proteomics_data <- proteomics_data %>%
      mutate(imputeVal = rowMeans(., na.rm = TRUE)) %>%
      mutate(., across(!imputeVal, ~ coalesce(.x, imputeVal))) %>%
      within(., rm('imputeVal'))
  } else if (imputation_method == "median") {
    # analogous to `imputation_method == "median"`. Uses `dplyr::rowwise` with
    # `dplyr::c_across()` in the `median()` call to achieve the same thing,
    # except calculating the medians of a given feature across all samples.
    proteomics_data <- proteomics_data %>%
      rowwise %>%
      mutate(imputeVal = median(c_across(), na.rm = TRUE)) %>%
      mutate(., across(!imputeVal, ~ coalesce(.x, imputeVal))) %>%
      within(., rm('imputeVal'))
  } else {
    stop()
  }

  # finally we store the imputed data back into the `SummarizedExperiment`
  # object and return it.
  SummarizedExperiment::assays(
    experiment,
    withDimnames = FALSE
  )[[assay]] <- proteomics_data
  return(experiment)
}


filter_for_feature_missingness <- function(experiment, missingness_cutoff = 0) {
  pdat <- assay(experiment)

  # TODO: add check to make sure missgness_cutoff is [0, 1]
  pdat %<>%
    .[rowSums(!is.na(.)) / ncol(.) >= missingness_cutoff, ]

  if (ncol(pdat) == 0) {
    print('All data removed')
    return(NULL)
  }
  rdat <- rowData(experiment)[rownames(pdat), ]

  experiment <- SummarizedExperiment(
    assays = list(proteomics = pdat),
    rowData = rdat,
    colData = colData(experiment)
  )
  return(experiment)
}


filter_for_sample_cutoff <- function(experiment, samp_cutoff = 0.1) {
  pdat <- assay(experiment)

  #What if we add a sample cutoff?
  pdat %<>%
    .[, colSums(!is.na(.)) / nrow(.) >= samp_cutoff] # %>%
  #rownames_to_column(., var = "Protein.Group")

  rdat <- colData(experiment)[colnames(pdat), ]

  experiment <- SummarizedExperiment(
    assays = list(proteomics = pdat),
    rowData = rowData(experiment),
    colData = rdat
  )

  return(experiment)
}
