#' Proteomics preprocessing utilities
#'
#' This file contains helper functions for preparing raw proteomics input data
#' into a standardized `SummarizedExperiment` object suitable for downstream
#' filtering, imputation, transformation, and statistical analysis. The
#' functions cover metadata cleaning, global and phosphoproteomics data
#' formatting, missingness filtering, imputation, and log transformation.
#'
#' @description
#' The file is organized around a small preprocessing pipeline for proteomics
#' data. It cleans sample metadata, converts raw abundance tables into a common
#' feature-by-sample structure, aligns annotations with assay data, optionally
#' removes low-quality features or samples, imputes missing values, and log2-
#' transforms the final assay data before downstream modeling.
#'
#' @details
#' Core functions in this file include `prepare_proteomics_data()`,
#' `clean_meta_data()`, `clean_global_proteomics_data()`,
#' `clean_phospho_proteomics_data()`, `log_transform()`,
#' `impute_proteomics()`, `filter_for_feature_missingness()`, and
#' `filter_for_sample_cutoff()`. Together they standardize proteomics inputs
#' for analysis in a reproducible and package-friendly workflow.
#'
#' @keywords internal
#' @name proteomics_data_prep

library(tibble)
library(dplyr)
library(tidyr)
library(magrittr)

library(SummarizedExperiment)
library(limma)
library(matrixStats)
library(sva)


#' Prepare proteomics data into a SummarizedExperiment
#'
#' Cleans and normalizes proteomics input tables, aligns metadata, and returns
#' a `SummarizedExperiment` ready for downstream filtering, imputation, and
#' statistical analysis. The function supports both global and phosphoproteomics
#' workflows and applies optional sample/feature filtering and log
#' transformation.
#'
#' @param proteomics_data A data frame or tibble containing proteomics intensity
#'   values, typically with a feature identifier column such as protein or
#'   phosphosite names.
#' @param proteomics_meta A data frame containing sample metadata. The metadata
#'   must include a `condition` field that can be parsed into group and
#'   replicate information.
#' @param type Character scalar indicating the assay type. Must be either
#'   `"global"` or `"phospho"`.
#' @param ... Additional optional arguments controlling preprocessing behavior:
#'   \itemize{
#'     \item{`sample_col_prefix`}{Character string used to identify the sample
#'       columns in the raw proteomics table. Defaults to `"cNF"`.}
#'     \item{`sample_missingness_cutoff`}{Numeric value in `[0, 1]` specifying
#'       the minimum fraction of non-missing values required for each sample to
#'       be retained. Mutually exclusive with `samples_to_exclude`.}
#'     \item{`samples_to_exclude`}{Character vector of sample names to remove
#'       prior to downstream analysis. Mutually exclusive with
#'       `sample_missingness_cutoff`.}
#'     \item{`feature_missingness_cutoff`}{Numeric value in `[0, 1]` specifying
#'       the minimum fraction of non-missing values required for each feature to
#'       be retained.}
#'     \item{`imputation_method`}{Character string naming the imputation method
#'       to apply before log transformation. Supported values are `"SampMin"`,
#'       `"means"`, and `"median"`.}
#'     \item{`group_contains_spheroid`}{Logical scalar indicating whether the
#'       metadata `condition` values include a spheroid identifier, which changes
#'       how the group field is parsed during metadata cleaning.}
#'   }
#'
#' @return A `SummarizedExperiment` object with assay values stored in
#'   `assay(experiment, "proteomics")`, feature metadata in `rowData()`, and
#'   sample metadata in `colData()`.
#'
#' @details The function validates the selected assay type, cleans both the data
#' and metadata, constructs the `SummarizedExperiment`, applies optional missing-
#' value filtering, imputes missing values when requested, and then log2-
#' transforms the expression matrix. Errors are raised when incompatible
#' filtering arguments are supplied or when required metadata is missing.
#'
#' @examples
#' \dontrun{
#' se <- prepare_proteomics_data(
#'   proteomics_data = protein_data,
#'   proteomics_meta = sample_metadata,
#'   type = "global",
#'   sample_col_prefix = "cNF",
#'   sample_missingness_cutoff = 0.5,
#'   feature_missingness_cutoff = 0.2,
#'   imputation_method = "median",
#'   group_contains_spheroid = FALSE
#' )
#'
#' assay(se)
#' }
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

#' Clean and parse metadata for proteomics sample grouping
#'
#' Splits the `condition` field in sample metadata into standardized group and
#' replicate columns, and optionally separates spheroid information when the
#' condition labels include it. The cleaned metadata is returned with sample
#' names stored as row names for use in a `SummarizedExperiment`.
#'
#' @param data A data frame containing sample metadata, including a `condition`
#'   column.
#' @param group_contains_spheroid Logical scalar indicating whether the condition
#'   string includes a spheroid identifier before the timepoint and agent fields.
#'   Defaults to `FALSE`.
#'
#' @return A cleaned metadata tibble or data frame with sample names as row
#'   names and parsed columns such as `group`, `replicate`, `timepoint`,
#'   `agent`, and optionally `spheroid`.
#'
#' @details The function uses `separate_wider_regex()` to split the `condition`
#' field into a `group` and `replicate` component, then further separates the
#' group value into either `timepoint` and `agent` or `spheroid`, `timepoint`,
#' and `agent` depending on the value of `group_contains_spheroid`.
#'
#' @examples
#' \dontrun{
#' metadata <- data.frame(
#'   sample = c("sample_1", "sample_2"),
#'   condition = c("control_1", "treated_2")
#' )
#'
#' clean_meta_data(metadata)
#' }
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

#' Clean and format global proteomics data
#'
#' Extracts sample intensity columns from a raw global proteomics table,
#' keeps the feature metadata needed for downstream analysis, and returns a
#' cleaned data object together with a mapping table for feature annotations.
#'
#' @param data A data frame containing raw global proteomics measurements.
#'   The table must include a `Protein.Group` identifier column and one or more
#'   sample columns whose names contain the prefix specified by
#'   `sample_col_prefix`.
#' @param sample_col_prefix Character string used to select the sample columns
#'   from the raw data. Defaults to `"cNF"`.
#'
#' @return A named list with two elements:
#' \itemize{
#'   \item{`data`}{A cleaned tibble/data frame with `Protein.Group` set as the
#'     row name and only the selected sample columns retained.}
#'   \item{`mapping`}{A feature annotation table describing the retained
#'     proteins, including parsed gene information.}
#' }
#'
#' @details This helper keeps the expression matrix aligned with the protein
#' metadata by selecting only the sample columns matching the sample prefix and
#' then returning the corresponding feature annotation table. The `Genes` column
#' is split on semicolons so the mapping table retains a simplified annotation
#' for each protein.
#'
#' @examples
#' \dontrun{
#' prot_data <- clean_global_proteomics_data(
#'   data = global_protein_table,
#'   sample_col_prefix = "cNF"
#' )
#'
#' prot_data$data[1:5, 1:3]
#' }
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

#' Clean and format phosphoproteomics data
#'
#' Builds a phosphosite-level expression matrix from a raw phosphoproteomics
#' table, creates a unique feature identifier for each phosphosite, and returns
#' both the cleaned abundance matrix and a mapping table with annotation data.
#'
#' @param data A data frame containing phosphoproteomics measurements. It must
#'   include columns such as `Protein`, `Residue`, `Site`, and `Gene.Names`,
#'   along with sample intensity columns matching `sample_col_prefix`.
#' @param sample_col_prefix Character string used to select the sample columns
#'   from the raw data. Defaults to `"cNF"`.
#'
#' @return A named list with two elements:
#' \itemize{
#'   \item{`data`}{A cleaned tibble/data frame with phosphosite identifiers set
#'     as row names and only the selected sample columns retained.}
#'   \item{`mapping`}{A feature annotation table containing phosphosite
#'     metadata, gene annotation, median abundance, and a quantile rank.}
#' }
#'
#' @details The function combines protein and site information into a unique
#' `Protein.With.Phospho.Site` identifier, removes redundant annotation columns,
#' and computes median abundance and quantile ranking for each phosphosite. This
#' keeps the expression matrix aligned with the phosphosite annotation table.
#'
#' @examples
#' \dontrun{
#' phospho_data <- clean_phospho_proteomics_data(
#'   data = phospho_table,
#'   sample_col_prefix = "cNF"
#' )
#'
#' phospho_data$data[1:5, 1:3]
#' }
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


#' Log-transform an assay in a SummarizedExperiment
#'
#' Applies a base-2 log transformation to the selected assay within a
#' `SummarizedExperiment`, converting abundance values to log-scale for
#' downstream statistical modeling and visualization.
#'
#' @param experiment A `SummarizedExperiment` object containing the assay to
#'   transform.
#' @param assay Character scalar naming the assay to transform. Defaults to
#'   `"proteomics"`.
#'
#' @return The original `SummarizedExperiment` object with the specified assay
#'   replaced by its log2-transformed values.
#'
#' @details The function coerces all values in the selected assay to numeric and
#' then applies `log2()`. The assay name is preserved so downstream code can
#' continue to reference the same assay slot.
#'
#' @examples
#' \dontrun{
#' se <- log_transform(se, assay = "proteomics")
#' assay(se)[1:5, 1:3]
#' }
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

#' Impute missing proteomics values in a SummarizedExperiment assay
#'
#' Replaces missing values in a selected assay using one of several supported
#' imputation strategies, then returns the updated `SummarizedExperiment` object.
#' This is useful when preparing proteomics data for downstream analyses that do
#' not tolerate missing values.
#'
#' @param experiment A `SummarizedExperiment` object containing the assay to
#'   impute.
#' @param assay Character scalar naming the assay to impute. Defaults to
#'   `"proteomics"`.
#' @param imputation_method Character string specifying the imputation approach.
#'   Supported values are `"SampMin"`, `"means"`, and `"median"`.
#'
#' @return The original `SummarizedExperiment` object with missing values
#'   replaced according to the selected imputation method.
#'
#' @details The function coerces assay values to numeric and then fills missing
#' entries according to the selected strategy: `SampMin` uses each sample's
#' minimum observed value, `means` uses the mean intensity for each feature, and
#' `median` uses the median intensity for each feature.
#'
#' @examples
#' \dontrun{
#' se <- impute_proteomics(
#'   experiment = se,
#'   assay = "proteomics",
#'   imputation_method = "median"
#' )
#' }
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

#' Filter features by missingness threshold
#'
#' Removes features whose proportion of non-missing values falls below a user-
#' specified cutoff. This is useful for excluding low-quality proteins or
#' phosphosites before downstream statistical analysis.
#'
#' @param experiment A `SummarizedExperiment` object containing the assay to
#'   filter.
#' @param missingness_cutoff Numeric value in `[0, 1]` representing the minimum
#'   fraction of non-missing values required for each feature to be retained.
#'   Defaults to `0`.
#'
#' @return A filtered `SummarizedExperiment` object containing only features that
#'   pass the missingness threshold. If no features remain, `NULL` is returned.
#'
#' @details The function computes the proportion of non-missing values for each
#'   row and retains only features with a proportion greater than or equal to
#'   `missingness_cutoff`. The row metadata is subset to match the retained
#'   features.
#'
#' @examples
#' \dontrun{
#' se_filtered <- filter_for_feature_missingness(se, missingness_cutoff = 0.7)
#' }
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

#' Filter samples by missingness threshold
#'
#' Removes samples whose proportion of non-missing values falls below a user-
#' specified cutoff. This helps exclude low-quality or poorly measured samples
#' before downstream differential expression and visualization steps.
#'
#' @param experiment A `SummarizedExperiment` object containing the assay to
#'   filter.
#' @param samp_cutoff Numeric value in `[0, 1]` representing the minimum
#'   fraction of non-missing values required for each sample to be retained.
#'   Defaults to `0.1`.
#'
#' @return A filtered `SummarizedExperiment` object containing only samples that
#'   pass the missingness threshold.
#'
#' @details The function computes the proportion of non-missing values for each
#'   sample and discards those below the cutoff. The sample metadata is subset to
#'   match the retained samples.
#'
#' @examples
#' \dontrun{
#' se_filtered <- filter_for_sample_cutoff(se, samp_cutoff = 0.5)
#' }
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
