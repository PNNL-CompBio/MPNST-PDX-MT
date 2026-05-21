
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
    as_tibble(.)

  # This makes sure that the differential expression results for each feature
  # are properly aligned

  SummarizedExperiment::rowData(experiment) <- left_join(
      as_tibble(SummarizedExperiment::rowData(experiment)),
      result_table,
      by = align_by)

  return(experiment)
}

prep_prot_data <- function(proteomics_data, proteomics_meta,
                           type = c('global', 'phospho'),
                           missingness_cutoff = 0.5,
                           samp_cutoff = 0.9,
                           sample_col_prefix = "cNF"){

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
    ) %>%
    separate_wider_delim(
      .,
      col = "group",
      delim = "_",
      names = c("timepoint", "agent"),
      too_few = "align_start",
      cols_remove = FALSE
    ) |>
    tibble::column_to_rownames('sample') # adding rownames for SummarizedExperiment

  if(sel_type == 'global'){
    proteomics_data_clean <- proteomics_data %>%
      column_to_rownames(., var = "Protein.Group") %>%
      select(., contains(sample_col_prefix))

    #What if we add a sample cutoff?
    proteomics_data_clean %<>%
      .[,colSums(!is.na(.))/nrow(.) >= samp_cutoff]# %>%
      #rownames_to_column(., var = "Protein.Group")

    # TODO: add check to make sure missgness_cutoff is [0, 1]
    proteomics_data_clean %<>%
      .[rowSums(!is.na(.))/ncol(.) >= missingness_cutoff,] %>%
      rownames_to_column(., var = "Protein.Group")

    proteomicsMapping <- proteomics_data[0:4] %>%
      separate_wider_delim(
        .,
        col = "Genes",
        names = c("Genes", NA),
        delim = ";",
        cols_remove=TRUE,
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
      within(., rm("Protein.Names", "Gene.Names", "Sequence"))

    ##add in missingness here
    proteomics_data_clean %<>%
      .[,colSums(!is.na(.))/nrow(.) >= samp_cutoff]#

    # TODO: add check to make sure missgness_cutoff is [0, 1]
    proteomics_data_clean %<>%
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


#-------------------
# plotting functions
#-------------------

color_vals_grey_first <- c('#BBBBBB', '#EE7733', '#0077BB', '#33BBEE', '#EE3377', '#CC3311', '#009988')
color_vals_black_first <- c('#000000', '#EE7733', '#0077BB', '#33BBEE', '#EE3377', '#CC3311', '#009988')
color_vals_default <- c('#EE7733', '#0077BB', '#33BBEE', '#EE3377', '#CC3311', '#009988', '#BBBBBB')


plot_pca <- function(experiment, agent, type = c('global', 'phospho'), missingness = 'None', imputation = "None"){

  prot_type <- match.arg(type)

  plates <- unique(colData(experiment[, experiment$agent ==  agent])$plate)
  plates <- append(plates, 10)
  ctrl <- ifelse(test = (agent == 'Trab'), yes = "Water", no = "DMSO")
  sel <- c(ctrl, agent)



  experiment_subset <- experiment[, (experiment$plate %in% plates) & experiment$agent %in% sel]

  principal_components <- prcomp(t(assay(experiment_subset))) %$%
    cbind(x, colData(experiment_subset))

  plot_title <- base::paste("Spheroid", prot_type, "proteomics for:", agent)
  plot_subtitle <- base::paste("Missingness filter:", missingness ,"; Imputation:", imputation)
  color_vals <- c("#66c2a5", "#fc8d62", "#8da0cb")
  colors <- setNames(color_vals, c(agent, "DMSO", "Water"))

  plot <- (
    ggplot(
      principal_components,
      aes(
        x = PC1,
        y = PC2,
        # col = timepoint,
        # shape = agent,
        col = agent,
        shape = timepoint,
        label = condition
      )
    )
    + geom_point(size = 2)
    + geom_text_repel(max.overlaps = 10, size = 2)
    + theme_minimal()
    + theme(legend.position = "bottom")
    + scale_color_manual(values = colors)
    # + scale_color_manual(values = c("DMSO" = "#e66101", "Mirda" = "#5e3c99", "NA" = "grey"))
    + labs(
      title = plot_title,
      subtitle = "50% missingness filter; SampMin imputation"
    )

  )
  return(plot)


}


plot_pca_by_plate <- function(experiment, plate, type = c('global', 'phospho'), missingness = 'None', imputation = "None"){

  prot_type <- match.arg(type)

  experiment_subset <- experiment[, experiment$plate == plate]

  principal_components <- prcomp(t(assay(experiment_subset))) %$%
    cbind(x, colData(experiment_subset))

  plot_title <- base::paste("Spheroid", prot_type, "proteomics for Plate:", plate)
  plot_subtitle <- base::paste("Missingness filter:", missingness ,"; Imputation:", imputation)
  # color_vals <- c("#66c2a5", "#fc8d62", "#8da0cb")

  agents <- unique(colData(experiment_subset)$agent)
  ctrl <- ifelse(test = ("DMSO" %in% agents), yes = "DMSO", no = "Water")
  agents_sel <- agents[agents != ctrl]
  colors <- setNames(color_vals_black_first, append(c(ctrl), agents_sel))

  plot <- (
    ggplot(
      principal_components,
      aes(
        x = PC1,
        y = PC2,
        # col = timepoint,
        # shape = agent,
        col = agent,
        shape = timepoint,
        label = condition
      )
    )
    + geom_point(size = 2)
    + geom_text_repel(max.overlaps = 10, size = 2)
    + theme_minimal()
    + theme(legend.position = "bottom")
    + scale_color_manual(values = colors)
    # + scale_color_manual(values = c("DMSO" = "#e66101", "Mirda" = "#5e3c99", "NA" = "grey"))
    + labs(
      title = plot_title,
      subtitle = plot_subtitle
    )

  )
  return(plot)


}


is_outlier <- function(x) {
  return(
    x < quantile(x, 0.25) - 1.5 * IQR(x) |
    x > quantile(x, 0.75) + 1.5 * IQR(x)
  )
}

plot_per_plate_na_overview <- function(experiment, type = c('global', 'phospho')){

  prot_type <- match.arg(type)

  summarized_data <- as_tibble(assay(experiment)) %>%
    summarise(across(everything(), ~ sum(is.na(.)))) %>%
    pivot_longer(cols = everything(), names_to = "sample_id", values_to = "NAs")

  metadata <- colData(experiment) %>%
    as.data.frame(.) %>%
    rownames_to_column(., var = "sample_id") %>% as_tibble(.)

  data_plot <- inner_join(metadata, summarized_data, by='sample_id') %>%
    group_by(plate) %>%
    mutate(outlier = ifelse(is_outlier(NAs), condition, NA)) %>%
    mutate(., agent = replace_when(agent,
                                   agent %in% c("DMSO", "Water") ~ "DMSO/Water"))


  plot_title <- base::paste("Outliers in spheroid", prot_type, "proteomics")
  plot_subtitle <- base::paste("# of NAs per sample out of", dim(assay(experiment))[1], ifelse(prot_type == "global", "proteins", "phospho-sites"))

  plot <- (
    ggplot(data_plot, aes(
      x = plate,
      y = NAs,
      group = plate,
      xmin = 1,
      xmax = 10
    ))
    + geom_boxplot()
    + geom_label_repel(aes(label = outlier),
                       na.rm = TRUE,
                       max.overlaps = 10,
                       size = 2,
                       min.segment.length = 0)
    + theme_minimal()
    + scale_x_continuous(breaks = seq(1, 10, 1))
    + labs(
      title = plot_title,
      subtitle = plot_subtitle
    )
  )

  return(plot)
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


plot_diff_ex <- function(experiment, agent, timepoint, log_fc = 1, pval = 0.05){

  ctrl <- ifelse(test = (agent == 'Trab'), yes = "Water", no = "DMSO")
  group_1 = paste0(timepoint, "_", ctrl)
  group_2 = paste0(timepoint, "_", agent)

  plate_sel <- unique(experiment[, experiment$group == group_2]$plate)
  experiment_subset <- experiment[, experiment$plate == plate_sel]

  align_by_group <- ifelse(test = (prot_type == 'global'), yes = "Protein.Group", no = "Protein.With.Phospho.Site")
  diff_ex <- calc_diff_ex(experiment_subset, group_1, group_2, align_by=align_by_group)

  if(prot_type == "global"){
    dfPlot <- extract_diff_ex_features(diff_ex, "Genes", log_fc, pval)
  } else if(prot_type == "phospho"){
    dfPlot <- extract_diff_ex_features(diff_ex, "Gene.With.Phospho.Site", log_fc, pval)
  } else {
    stop(paste0("Invalid prot_type selected: ", prot_type))
  }

  # plot_title = paste('Spheroid', prot_type ,'proteomics:', agent, "vs.", ctrl, "at", timepoint)
  plot_subtitle = base::paste0(
    timepoint, "; ",
    "#down: ", dim(dfPlot[dfPlot$expression_change == 'down',])[1], " / ",
    "#up: ", dim(dfPlot[dfPlot$expression_change == 'up',])[1]
  )

  plot <- (
    ggplot(
      data = dfPlot,
      aes(
        x = logFC,
        y = -log10(adj.P.Val),
        col = expression_change,
        label = labels
      )
    )
    + geom_vline(xintercept = c(-log_fc, log_fc), col = "gray", linetype = 'dashed')
    + geom_hline(yintercept = -log10(pval), col = "gray", linetype = 'dashed')
    + geom_point(size = .75)
    + geom_text_repel(max.overlaps = 5, size = 2)
    + scale_color_manual(values = c('up' = '#74add1', 'down' = '#f46d43', 'not significant' = 'grey'))
    + theme_minimal()
    + theme(legend.position = "bottom")
    + labs(
      # title = plot_title,
      subtitle = plot_subtitle,
      col = 'expression',
      y = expression("-log"[10]*"p-Value")
    )
  )

  plot
}

plot_diff_ex_grid <- function(plots, plot_title=NULL){

  plot_grid <- ggpubr::ggarrange(plotlist = plots, ncol = 3, common.legend = TRUE, legend = "bottom")
  plot_annotated <- ggpubr::annotate_figure(plot_grid, top = ggpubr::text_grob(plot_title, size = 14))

  plot_annotated
}
