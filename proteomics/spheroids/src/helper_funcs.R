library(tibble)
library(dplyr)
library(tidyr)
library(magrittr)

library(SummarizedExperiment)
library(limma)
library(matrixStats)
library(sva)


# DEPRECATED!
# Currently in here for backwards compatibility, with MN2 & JH097c notebooks
# TODO: remove after fixing notebooks
prep_experiment_data <- function(
  prot_data,
  meta_data,
  missingness_cutoff = 0,
  samp_cutoff = list(global = 0, phospho = 0),
  imputation_method = c("none", 'SampMin', 'means', 'median'),
  batch_correct = FALSE,
  group_contains_spheroid = FALSE
) {
  imputation_method_sel <- match.arg(imputation_method)

  if (
    !is.na(missingness_cutoff) &
      !(missingness_cutoff >= 0 | missingness_cutoff <= 1)
  ) {
    stop(
      base::paste(
        "missingness_cutoff needs to be '0 <= missingness_cutoff <=1'",
        "or 'NULL'. Supplied value:",
        missingness_cutoff
      ),
      call. = FALSE
    )
  }

  experiments <- list()
  # for(data_type in c('global', 'phospho')){
  for (data_type in c('global')) {
    experiment <- prep_prot_data(
      proteomics_data = prot_data[[data_type]],
      proteomics_meta = meta_data[[data_type]],
      missingness_cutoff = missingness_cutoff,
      samp_cutoff = samp_cutoff[[data_type]],
      sample_col_prefix = 'G_',
      type = data_type,
      imputation_method = imputation_method_sel,
      group_contains_spheroid = group_contains_spheroid
    )

    if (imputation_method_sel == "none") {
      experiment %<>%
        log_transform(., assay = "proteomics")
    }

    if (!(imputation_method_sel == 'none')) {
      experiment %<>%
        impute_proteomics(., assay = "proteomics", imputation_method_sel) %>%
        log_transform(., assay = "proteomics")
    }

    if (batch_correct && imputation_method_sel != 'none') {
      #first we have to remove zero hour beacuse it is confounded by batch
      #experiment <- experiment[,-which(colData(experiment)$timepoint=='0hr')]
      prot_matrix <- assay(experiment)
      batch_var <- colData(experiment)$plate
      cond_var <- colData(experiment)$agent
      mod_matrix <- model.matrix(~ as.factor(agent), data = colData(experiment))
      corrected_prot <- sva::ComBat(
        dat = prot_matrix,
        batch = batch_var,
        mod = mod_matrix,
        par.prior = FALSE
      )
      assay(experiment, withDimnames = FALSE) <- corrected_prot
    }
    experiments[[data_type]] <- experiment
  }

  experiments
}


# DEPRECATED!
# Currently in here for backwards compatibility, with MN2 & JH097c notebooks
# TODO: remove after fixing notebooks
prep_prot_data <- function(
  proteomics_data,
  proteomics_meta,
  type = c('global', 'phospho'),
  missingness_cutoff = 0.0,
  samp_cutoff = 0,
  sample_col_prefix = "cNF",
  imputation_method = 'none',
  group_contains_spheroid = FALSE
) {
  # checking that either global or phospho is selected for type
  # match.arg() will throw an error otherwise and execution will stop
  sel_type <- match.arg(type)

  # cleaning up the metadata by
  #
  #
  if (group_contains_spheroid == TRUE) {
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
        names = c("spheroid", "timepoint", "agent"),
        too_few = "align_start",
        cols_remove = FALSE
      ) |>
      tibble::column_to_rownames('sample')
  } else {
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
  }

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
  } else if (sel_type == 'phospho') {
    proteomics_data_clean <- proteomics_data %>%
      unite(., "Phospho.Site", Residue, Site, sep = "", remove = TRUE) %>%
      unite(
        .,
        "Protein.With.Phospho.Site",
        Protein,
        Phospho.Site,
        sep = "-"
      ) %>%
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
      semi_join(., proteomics_data_clean, by = "Protein.With.Phospho.Site")

    proteomics_data_clean %<>%
      column_to_rownames(., var = "Protein.With.Phospho.Site")

    #calculate median abundnace for each phosphosite and rank
    medAbund <- apply(proteomics_data_clean, 1, function(x) {
      median(x, na.rm = T)
    })
    qvals <- quantile(medAbund, probs = c(0.1, 0.2, 0.5, 1), na.rm = TRUE)
    siteQuantile <- sapply(medAbund, function(x) {
      names(qvals)[which(x <= qvals)[1]]
    })

    proteomicsMapping$medAbundance <- medAbund
    proteomicsMapping$quantile <- siteQuantile
  }
  experiment <- SummarizedExperiment(
    assays = list(proteomics = proteomics_data_clean),
    rowData = proteomicsMapping,
    colData = proteomics_meta_clean[colnames(proteomics_data_clean), ]
  )

  experiment <- filter_for_sample_cutoff(experiment, samp_cutoff)
  experiment <- filter_for_feature_missingness(experiment, missingness_cutoff)
  return(experiment)
}

#-------------------
# stats functions
#-------------------

exp_subset <- function(experiment, agent, timepoint) {
  ctrl <- ifelse(test = (agent == 'Trab'), yes = "Water", no = "DMSO")
  group_1 <- paste0(timepoint, "_", ctrl)
  group_2 <- paste0(timepoint, "_", agent)
  ##the old data doesn't have a plate
  if ('plate' %in% names(colData(experiment))) {
    plate_sel <- unique(experiment[, experiment$group == group_2]$plate)
    experiment_subset <- experiment[, experiment$plate == plate_sel]
  } else {
    experiment_subset <- experiment[, experiment$timepoint == timepoint]
  }
  #experiment_subset <- experiment[, experiment$timepoint == timepoint]
  experiment_subset <- experiment_subset[,
    experiment_subset$agent %in% c(agent, ctrl)
  ]

  ##make sure that there are no additional features that we are adding back
  all_nas <- which(apply(assay(experiment_subset), 1, function(x) {
    all(is.na(x))
  }))
  #print(paste('Found',length(all_nas),'rows that are all na'))
  if (length(all_nas) > 0) {
    experiment_subset <- experiment_subset[-all_nas, ]
  }
  return(experiment_subset)
}

##plot enrichment results of a single drug across time points
format_enrich_result <- function(
  enrichlist, #list of enrichment results
  timepoints, #list of timepoints
  agent,
  pval_thresh = 0.05
) {
  ##get a data frame of all kinases
  names(enrichlist) <- timepoints

  all.kin <- do.call(
    rbind,
    lapply(timepoints, function(tp) {
      lres <- enrichlist[[tp]]
      data.frame(
        feature = lres$feature,
        zscore = lres$zscore,
        pval = lres$BH_pvalue,
        timepoint = tp
      )
    })
  )

  if (nrow(all.kin) == 0) {
    print(paste("No enrichmed features for", agent))
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
