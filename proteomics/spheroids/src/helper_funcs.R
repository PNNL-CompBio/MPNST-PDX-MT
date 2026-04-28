
log_transform <- function(experiment,
                          assay = "proteomics"){

  proteomics_data <- SummarizedExperiment::assays(experiment)[[assay]] %>%
    mutate(across(everything(), as.numeric)) %>% # turn everything numeric
    log2(.) # log transform all cells

  SummarizedExperiment::assays(
    experiment,
    withDimnames = FALSE)[[assay]] <- proteomics_data
  return(experiment)
}

impute_proteomics <- function(experiment,
                              assay = 'proteomics',
                              imputation_method = 'SampMin'){

  # retrieve assay from SummarizedExperiment and make sure that all values
  # are numeric such that we can impute
  proteomics_data <- SummarizedExperiment::assays(experiment)[[assay]] %>%
    mutate(across(everything(), as.numeric)) # turn every cell numeric

  if(imputation_method == "SampMin") {

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
      mutate(., across(!imputeVal, ~coalesce(.x, imputeVal))) %>%
      within(., rm('imputeVal'))

  } else if (imputation_method == "median") {

    # analogous to `imputation_method == "median"`. Uses `dplyr::rowwise` with
    # `dplyr::c_across()` in the `median()` call to achieve the same thing,
    # except calculating the medians of a given feature across all samples.
    proteomics_data <- proteomics_data %>%
      rowwise %>%
      mutate(imputeVal = median(c_across(), na.rm = TRUE)) %>%
      mutate(., across(!imputeVal, ~coalesce(.x, imputeVal))) %>%
      within(., rm('imputeVal'))

  } else {
    stop()
  }

  # finally we store the imputed data back into the `SummarizedExperiment`
  # object and return it.
  SummarizedExperiment::assays(
    experiment,
    withDimnames = FALSE)[[assay]] <- proteomics_data
  return(experiment)

}

calc_diff_ex <- function(experiment,
                         group_1,
                         group_2,
                         assay = "proteomics",
                         align_by = "Protein.Group"){


  # making sure that the 'align_by' group is present in the SummarizedExperiment
  if(!align_by %in% colnames(rowData(experiment))){
    stop(
      base::paste(align_by,
             "not in rowData. Available columns:",
             paste0(colnames(rowData(experiment)), collapse = ","),
             sep = " "
             ),
      call. = FALSE)
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
  dat <- SummarizedExperiment::assays(experiment)[[assay]]

  result_table <- dat[, c(samples_2, samples_1)] %>%
    limma::lmFit(., design) %>%
    limma::eBayes(.) %>%
    limma::topTable(., coef = 2, number = Inf, sort.by = "none")

  # The assumption is that the input `SummarizedExperiment` had only one assay
  # associated with it, hence we will just augment the object and return it
  # containing the results from limma.

  result_table <- result_table %>%
    rownames_to_column(., var=align_by) %>%
    as.tibble(.)

  # This makes sure that the differential expression results for each feature
  # are properly aligned

  SummarizedExperiment::rowData(experiment) <- left_join(
      as.tibble(SummarizedExperiment::rowData(experiment)),
      result_table,
      by = align_by)

  return(experiment)
}

prep_global_prot_data <- function(proteomics_data, proteomics_meta, missingness_cutoff = 0.5) {

  proteomicsMetaClean <- proteomics_meta %>%
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
    )|>
    tibble::column_to_rownames('sample') ##SG: rownames reflect sample identifier


  proteomicsDataClean <- proteomicsData %>%
    column_to_rownames(., var = "Protein.Group") %>%
    within(., rm("Protein.Names", "Genes", "First.Protein.Description"))

  # TODO: add check to make sure missgness_cutoff is [0, 1]
  proteomicsDataClean %<>%
    .[rowSums(!is.na(.))/ncol(.) >= missingness_cutoff,] %>%
    rownames_to_column(., var = "Protein.Group")

  proteomicsMapping <- proteomicsData[0:4] %>%
    separate_wider_delim(
      .,
      col = "Genes",
      names = c("Genes", NA),
      delim = ";",
      cols_remove=TRUE,
      too_few = "align_start",
      too_many = "drop"
    ) %>%
    semi_join(., proteomicsDataClean, by = "Protein.Group")

  proteomicsDataClean %<>% column_to_rownames(., var = "Protein.Group")

  experiment <- SummarizedExperiment(
    assays = list(proteomics = proteomicsDataClean),
    rowData = proteomicsMapping,
    colData = proteomicsMetaClean[colnames(proteomicsDataClean),]
  )

  return(experiment)
}

prep_phospho_prot_data <- function(proteomics_data, proteomics_meta, missingness_cutoff = 0.5) {

   proteomicsMetaClean <- proteomics_meta %>%
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
       )|>
     tibble::column_to_rownames('sample') ##SG: rownames reflect sample identifier

   proteomicsDataClean <- proteomics_data %>%
     unite(., "Phospho.Site", Residue, Site, sep="", remove=TRUE) %>%
     unite(., "Protein.With.Phospho.Site", Protein, Phospho.Site, sep = "-") %>%
     column_to_rownames(., var = "Protein.With.Phospho.Site") %>%
     within(., rm("Protein.Names", "Gene.Names", "Sequence"))

   # TODO: add check to make sure missgness_cutoff is [0, 1]
   proteomicsDataClean %<>%
     .[rowSums(!is.na(.))/ncol(.) >= missingness_cutoff,] %>%
     rownames_to_column(., var = "Protein.With.Phospho.Site")

   proteomicsMapping <- proteomics_data[0:5] %>%
     separate_wider_delim(
       .,
       col = "Gene.Names",
       names = c("Gene.Names", NA),
       delim = ";",
       cols_remove=TRUE,
       too_few = "align_start",
       too_many = "drop"
       ) %>%
     unite(., "Phospho.Site", Residue, Site, sep="", remove=TRUE) %>%
     unite(., "Protein.With.Phospho.Site", Protein, Phospho.Site,
           sep = "-", remove = FALSE) %>%
     unite(., "Gene.With.Phospho.Site", Gene.Names, Phospho.Site,
           sep = "-", remove = FALSE) %>%
     semi_join(., proteomicsDataClean, by = "Protein.With.Phospho.Site")

   proteomicsDataClean %<>% column_to_rownames(., var = "Protein.With.Phospho.Site")

   #calculate median abundnace for each phosphosite and rank
   medAbund <- apply(proteomicsDataClean, 1, function(x) median(x,na.rm = T))
   qvals <- quantile(medAbund,probs = c(0.1,0.2,0.5,1))
   siteQuantile <- sapply(medAbund, function(x) names(qvals)[which(x <= qvals)[1]])


   proteomicsMapping$medAbundance = medAbund
   proteomicsMapping$quantile = siteQuantile

   experiment <- SummarizedExperiment(
     assays = list(proteomics = proteomicsDataClean),
     rowData = proteomicsMapping,
     colData = proteomicsMetaClean[colnames(proteomicsDataClean),]
   )

   return(experiment)
}
