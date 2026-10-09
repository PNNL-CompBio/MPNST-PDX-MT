

library(rstatix)
log_transform <- function(experiment,
                          assay = "proteomics"){

  proteomics_data <- SummarizedExperiment::assay(experiment,assay) %>%
    mutate(across(everything(), as.numeric)) %>% # turn everything numeric
    log2(.) # log transform all cells

  SummarizedExperiment::assay(
    experiment,
    withDimnames = FALSE) <- proteomics_data
  SummarizedExperiment::assayNames(experiment) <- assay
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
                         assay_name = "proteomics",
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
  dat <- SummarizedExperiment::assays(experiment)[[assay_name]]

  result_table <- dat[, c(samples_2, samples_1)] %>%
    limma::lmFit(., design) %>%
    limma::eBayes(.) %>%
    limma::topTable(., coef = 2, number = Inf, sort.by = "none")

  # The assumption is that the input `SummarizedExperiment` had only one assay
  # associated with it, hence we will just augment the object and return it
  # containing the results from limma.

  result_table <- result_table %>%
    rownames_to_column(., var=align_by) %>%
    as_tibble(.)

  # This makes sure that the differential expression results for each feature
  # are properly aligned

  SummarizedExperiment::rowData(experiment) <- left_join(
      as_tibble(SummarizedExperiment::rowData(experiment)),
      result_table,
      by = align_by)

  return(experiment)
}

filter_for_feature_missingness <- function(experiment,
                                           missingness_cutoff = 0){

  pdat <- assay(experiment)

  # TODO: add check to make sure missgness_cutoff is [0, 1]
  pdat %<>%
    .[rowSums(!is.na(.))/ncol(.) >= missingness_cutoff,]

  if (ncol(pdat) == 0) {
    print('All data removed')
    return(NULL)
  }
  rdat <- rowData(experiment)[rownames(pdat),]

  experiment <- SummarizedExperiment(
    assays = list(proteomics = pdat),
    rowData = rdat,
    colData = colData(experiment)
  )
  return(experiment)
}


filter_for_sample_cutoff <- function(experiment,
                                     samp_cutoff=0.1){

  pdat <- assay(experiment)

  #What if we add a sample cutoff?
  pdat %<>%
    .[,colSums(!is.na(.))/nrow(.) >= samp_cutoff]# %>%
  #rownames_to_column(., var = "Protein.Group")

  rdat <- colData(experiment)[colnames(pdat),]

  experiment <- SummarizedExperiment(
    assays = list(proteomics = pdat),
    rowData = rowData(experiment),
    colData = rdat
  )

  return(experiment)
}


##this is a function to prep experiment data, calls the prep_prot_data below
##this is only applicable to the MN2 dataset
prep_experiment_data <- function(
    p_data,
    data_type = 'global',
    meta_data,
    missingness_cutoff=0,
    samp_cutoff = 0,
    imputation_method=c("none", 'SampMin', 'means', 'median'),
    batch_correct = FALSE
    ){

  imputation_method_sel <- match.arg(imputation_method)

  if(!is.na(missingness_cutoff) & !(missingness_cutoff >=0 | missingness_cutoff <= 1)){
    stop(
      base::paste(
        "missingness_cutoff needs to be '0 <= missingness_cutoff <=1'",
        "or 'NULL'. Supplied value:", missingness_cutoff),
      call. = FALSE
    )
  }

  experiments <- list()
  experiment <- prep_prot_data(proteomics_data = p_data,
                                 proteomics_meta = meta_data,
                                 missingness_cutoff = missingness_cutoff,
                                 samp_cutoff = samp_cutoff,
                                 type = data_type,
                                 imputation_method=imputation_method_sel)

    if(imputation_method_sel == "none"){
      experiment %<>%
        log_transform(., assay = "proteomics")
    }

    if(!(imputation_method_sel == 'none')){
      experiment %<>%
        impute_proteomics(., assay = "proteomics", imputation_method_sel) %>%
        log_transform(., assay = "proteomics")

    }

    if(batch_correct && imputation_method_sel != 'none'){
      #first we have to remove zero hour beacuse it is confounded by batch
      #experiment <- experiment[,-which(colData(experiment)$timepoint=='0hr')]
      prot_matrix <- assay(experiment)
      batch_var = colData(experiment)$plate
      cond_var = colData(experiment)$agent
      mod_matrix <- model.matrix(~ as.factor(agent), data = colData(experiment))
      corrected_prot <- sva::ComBat(
             dat = prot_matrix,
             batch = batch_var,
             mod = mod_matrix,
             par.prior = FALSE
              )
    assay(experiment,withDimnames = FALSE) <- corrected_prot


   }


  experiment
}


prep_prot_data <- function(proteomics_data, proteomics_meta,
                           type = c('global', 'phospho'),
                           missingness_cutoff = 0.0,
                           samp_cutoff = 0,
                           sample_col_prefix = "cNF",
                           imputation_method = 'none'){

  # checking that either global or phospho is selected for type
  # match.arg() will throw an error otherwise and execution will stop
  sel_type <- match.arg(type)

   # cleaning up the metadata by
  #
  proteomics_meta_clean <- proteomics_meta %>%
    separate_wider_regex(
      .,
      col = "condition",
      c("group" = ".*", "_", "replicate" = ".*"),
      cols_remove = FALSE
    )  %>%
    separate_wider_delim(
      .,
      col = "group",
      delim = "_",
      names = c("timepoint", "agent"),
      too_few = "align_start",
      cols_remove = FALSE
    ) |>
    tibble::column_to_rownames('sample') # adding rownames for SummarizedExperiment

  if (sel_type == 'global') {
    proteomics_data_clean <- proteomics_data %>%
      column_to_rownames(., var = "Protein.Group") %>%
      select(., contains(sample_col_prefix)) |>
      rownames_to_column("Protein.Group")


    proteomicsMapping <- proteomics_data[0:4] %>%
      separate_wider_delim(
        .,
        col = "Genes",
        names = c("Genes", NA),
        delim = ";",
        cols_remove = TRUE,
        too_few = "align_start",
        too_many = "drop"
      ) %>%
      semi_join(., proteomics_data_clean, by = "Protein.Group")

    proteomics_data_clean %<>% column_to_rownames(., var = "Protein.Group")

  } else if( sel_type == 'phospho'){
    proteomics_data_clean <- proteomics_data %>%
      unite(., "Phospho.Site", Residue, Site, sep="", remove=TRUE) %>%
      unite(., "Protein.With.Phospho.Site", Protein, Phospho.Site, sep = "-") %>%
      column_to_rownames(., var = "Protein.With.Phospho.Site") %>%
      within(., rm("Protein.Names", "Gene.Names", "Sequence")) |>
      rownames_to_column('Protein.With.Phospho.Site')

    # ##add in missingness here
    # proteomics_data_clean %<>%
    #   .[,colSums(!is.na(.))/nrow(.) >= samp_cutoff]#
    #
    # # TODO: add check to make sure missgness_cutoff is [0, 1]
    # proteomics_data_clean %<>%
    #   .[rowSums(!is.na(.))/ncol(.) >= missingness_cutoff,] %>%
    #   rownames_to_column(., var = "Protein.With.Phospho.Site")

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
      semi_join(., proteomics_data_clean, by = "Protein.With.Phospho.Site")

    proteomics_data_clean %<>% column_to_rownames(., var = "Protein.With.Phospho.Site")

    #calculate median abundnace for each phosphosite and rank
    medAbund <- apply(proteomics_data_clean, 1, function(x) median(x,na.rm = T))
    qvals <- quantile(medAbund,probs = c(0.1,0.2,0.5,1),na.rm=TRUE)
    siteQuantile <- sapply(medAbund, function(x) names(qvals)[which(x <= qvals)[1]])

    proteomicsMapping$medAbundance = medAbund
    proteomicsMapping$quantile = siteQuantile

  }
  experiment <- SummarizedExperiment(
    assays = list(proteomics = proteomics_data_clean),
    rowData = proteomicsMapping,
    colData = proteomics_meta_clean[colnames(proteomics_data_clean),]
  )

  experiment <- filter_for_sample_cutoff(experiment, samp_cutoff)
  experiment <- filter_for_feature_missingness(experiment, missingness_cutoff)
  return(experiment)

}

#-------------------
# stats functions
#-------------------


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
      explained_variance = as.data.frame(summary(.)$importance)[2,]
    )

  prc_res

}


extract_diff_ex_features <- function(diff_ex, feature_type, log_fc = 1, pval = 0.05){

  df_ret <- rowData(diff_ex) |> as_tibble()

  if(!feature_type %in% names(df_ret)){
    stop(paste0(
      "feature_type '",
      feature_type,
      "' not present in rowData of diff_ex",
      )
    )
  }

  df_ret %<>% dplyr::select(
      feature = feature_type,
      'logFC',
      'adj.P.Val'
    ) %>%
    mutate(expression_change = case_when(
      (logFC >= log_fc & adj.P.Val < pval) ~ 'up',
      (logFC <= -log_fc & adj.P.Val < pval) ~ 'down',
      .default = 'not significant',
      )
    ) %>%
    mutate(labels = case_when(
      expression_change %in% c('up', 'down') ~ feature,
      .default = NA
      )
    )

  df_ret
}

exp_subset <- function(experiment,
                       agent,
                       timepoint){

  ctrl <- ifelse(test = (agent == 'Trab'), yes = "Water", no = "DMSO")
  group_1 = paste0(timepoint, "_", ctrl)
  group_2 = paste0(timepoint, "_", agent)
  ##the old data doesn't have a plate
  if ('plate' %in% names(colData(experiment))) {
    plate_sel <- unique(experiment[, experiment$group == group_2]$plate)
    experiment_subset <- experiment[, experiment$plate == plate_sel]
  } else{
    experiment_subset <- experiment[, experiment$timepoint == timepoint]
  }
  #experiment_subset <- experiment[, experiment$timepoint == timepoint]
  experiment_subset <- experiment_subset[, experiment_subset$agent %in% c(agent, ctrl)]

  ##make sure that there are no additional features that we are adding back
  all_nas <- which(apply(assay(experiment_subset),1, function(x) all(is.na(x))))
  #print(paste('Found',length(all_nas),'rows that are all na'))
  if (length(all_nas) > 0){
    experiment_subset <- experiment_subset[-all_nas,]
  }
  return(experiment_subset)


}

##plot enrichment results of a single drug across time points
format_enrich_result <- function(enrichlist,  #list of enrichment results
                               timepoints, #list of timepoints
                               agent,
                               pval_thresh = 0.05){

  ##get a data frame of all kinases
  names(enrichlist) <- timepoints

  all.kin <- do.call (rbind, lapply(timepoints, function(tp) {
    lres <- enrichlist[[tp]]
     data.frame(feature = lres$feature, zscore = lres$zscore,
                      pval = lres$BH_pvalue, timepoint = tp)
  }))

  if (nrow(all.kin) == 0) {
    print(paste("No enrichmed features for",agent))
    return(NULL)
  }


  ##now only select those kinases that are significant in at least one time
  tps <- all.kin |>
    subset(pval < pval_thresh) |>
    arrange(rev(zscore))


  forder <- unique(tps$feature)
  torder <- timepoints

 all.kin <- all.kin |>
   subset(feature %in% tps$feature) |>
   mutate(sig = pval < pval_thresh) |>
   mutate(agent = agent)

 return(all.kin)
}
