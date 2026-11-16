## ---- sourcing-cross-notebook-functions ----

src_dir <- "../../src/"
source(paste0(src_dir, 'helper_funcs.R'))
source(paste0(src_dir, 'data_prep_funcs.R'))

## ---- library-loading ----

library(readxl)
library(tibble)
library(dplyr)
library(tidyr)
library(stringr)
library(SummarizedExperiment)


## ---- synapse-client-setup ----

library(reticulate)
py_require("synapseclient")
sc <- import("synapseclient")

## ---- synapse-data-import ----

# using the synapseclient python package to retrieve the data files from Synapse
syn <- sc$Synapse()
syn$login(silent = TRUE)

# getting the global proteomics abundance data
synapseFileID <- "syn77329971"
fileEntity <- syn$get(synapseFileID)
global_dat <- readr::read_tsv(
  fileEntity$path,
  na = c("0", "", NA), # the file stored on synapse has '0's instead of NA or ''
)

# getting the meta data
synapseFileID <- "syn77615964"
fileEntity <- syn$get(synapseFileID)
meta_dat <- read_xlsx(fileEntity$path, sheet = "Combined")

syn$logout()


## ---- proteomics-data-cleaning ----

# cleaning up the column names for the global proteomics data
names(global_dat) %<>%
  gsub("D:\\\\Data_D\\\\LeDay\\\\Spheroid_7863\\\\rawdata\\\\", "", .) %>%
  gsub("^[0-9]*_MPNST_spheroid_3Patient_", "", .) %>%
  gsub("_Maggie_50SPD_[0-9]+Aug[0-9]+_ES906-3408_[0-9]*.raw", "", .)


global_dat %<>%
  dplyr::select(
    -c('First.Protein.Description', 'N.Sequences', 'N.Proteotypic.Sequences')
  )


## ---- metadata-cleaning ----

# getting the new column names for the meta data
new_col_names <- c(
  global_sample_id = "Proteomics Sample ID",
  # phospho_sample_id = "Phosphoproteomics Sample ID",
  agent = "Drug",
  timepoint = "Time",
  plate = "Plate",
  spheroid = "Speroid"
)

# cleaning up sample ids in meta data
meta_dat$`Proteomics Sample ID` %<>% str_sub(., 23)

# renaming columns in meta data and creating new columns for downstream analysis
meta_dat %<>%
  dplyr::rename(., all_of(new_col_names)) %>%
  rowwise() %>%
  mutate(plate_by_spheroid = paste(spheroid, plate, sep = "_")) %>%
  mutate(agent = gsub(" ", "", agent), timepoint = paste0(timepoint, "hr")) %>%
  mutate(group = paste(spheroid, timepoint, agent, sep = "_")) %>%
  group_by(group) %>%
  mutate(replicate = paste0('r', 1:n())) %>%
  mutate(
    condition = paste(spheroid, timepoint, agent, replicate, sep = '_')
  ) %>%
  ungroup()

meta_data_global <- meta_dat %>%
  select(., global_sample_id, condition, plate_by_spheroid) %>%
  dplyr::rename(., sample = global_sample_id) %>%
  dplyr::rename(., plate = plate_by_spheroid)


## ---- optional-inclusion-of-mn2 ----

if (!eval_mn2) {
  cat("MN2 data will not be included in the analysis.")
} else {
  syn <- sc$Synapse()
  syn$login(silent = TRUE)

  # getting the global proteomics abundance data
  synapseFileID <- "syn74552028"
  fileEntity <- syn$get(synapseFileID)
  mn2_global_dat <- readr::read_tsv(
    fileEntity$path,
    na = c("0", "", NA), # the file stored on synapse has '0's instead of NA or ''
  )

  # getting the meta data
  synapseFileID <- "syn74755705"
  fileEntity <- syn$get(synapseFileID)
  mn2_meta_dat <- read_xlsx(fileEntity$path, sheet = 2)

  syn$logout()

  names(mn2_global_dat) <- names(mn2_global_dat) %>%
    gsub("D:\\\\Data_D\\\\LeDay\\\\spheroid_MN02\\\\rawdata\\\\", "", .) %>%
    gsub("_16Apr26_Monty_ES906_2014_50SPD[x]*.raw", "", .) %>%
    gsub("^[0-9]*_cNF_spheroid_", "G_", .)

  mn2_global_dat %<>%
    dplyr::select(
      -c('First.Protein.Description', 'N.Sequences', 'N.Proteotypic.Sequences')
    )

  new_col_names <- c(
    sample = 'Sample ID \r\n(abbrev)',
    global_sample_id = "Proteomics Sample ID",
    phospho_sample_id = "Phosphoproteomics Sample ID",
    agent = "Drug",
    timepoint = "Time",
    plate = "Plate"
  )

  mn2_meta_dat %<>%
    rowwise() %>%
    mutate(Drug = gsub(" ", "", Drug), Time = paste0(Time, "hr")) %>%
    mutate(spheroid = "MN2") %>%
    mutate(group = paste(Time, Drug, sep = "_")) %>%
    mutate(Plate = paste(spheroid, Plate, sep = "_")) %>%
    group_by(group) %>%
    mutate(replicate = paste0('r', 1:n())) %>%
    mutate(condition = paste(spheroid, Time, Drug, replicate, sep = '_')) %>%
    dplyr::rename(., all_of(new_col_names)) %>%
    ungroup()

  mn2_meta_dat <- mn2_meta_dat[
    mn2_meta_dat$plate %in% c("MN2_1", "MN2_4", "MN2_6"),
  ]

  mn2_meta_data_global <- mn2_meta_dat %>%
    select(., global_sample_id, condition, plate) %>%
    dplyr::rename(., sample = global_sample_id) %>%
    mutate(sample = gsub("cNF_spheroid_", "G_", sample))

  cols_to_keep <- c(
    'Protein.Group',
    'Protein.Names',
    'Genes',
    mn2_meta_data_global$sample
  )

  mn2_global_dat %<>%
    dplyr::select(all_of(cols_to_keep))

  global_dat %<>%
    full_join(
      .,
      mn2_global_dat,
      by = c('Protein.Group', 'Protein.Names', 'Genes')
    )

  meta_data_global %<>%
    rbind(
      .,
      mn2_meta_data_global
    )
}

## ---- creating-summarized-experiment ----

experiment <- prepare_proteomics_data(
  proteomics_data = global_dat,
  proteomics_meta = meta_data_global,
  type = 'global',
  sample_col_prefix = "G_",
  group_contains_spheroid = TRUE
)

## ---- calculate-missingness ----

pdat <- assay(experiment)
p_missingness <- enframe(
  rowSums(is.na(pdat)) / ncol(pdat),
  name = 'protein',
  value = 'missingness'
)

p_missingness$outlier <- get_zscore_outliers(
  p_missingness$missingness,
  sensitivity = 3,
  sided = 'right'
)

p_missingness$outlier_less_sensitive <- get_zscore_outliers(
  p_missingness$missingness,
  sensitivity = 2,
  sided = 'right'
)

cat(paste0(
  "Outliers (Z > 3): ",
  length(p_missingness$protein[p_missingness$outlier]),
  "\n",
  "Outliers (Z > 2): ",
  length(p_missingness$protein[p_missingness$outlier_less_sensitive])
))

z_score_based_missingness_cutoffs <- c(
  'Z > 3' = min(p_missingness$missingness[
    p_missingness$outlier
  ]),
  'Z > 2' = min(p_missingness$missingness[
    p_missingness$outlier_less_sensitive
  ])
)

cat(paste0(
  "Outliers (missingness > 75th percentile): ",
  dim(p_missingness[
    p_missingness$missingness > quantile(p_missingness$missingness)[4],
  ])[1]
))

## ---- missingness-by-sample ----

summarized_data <- as_tibble(assay(experiment)) %>%
  summarise(across(everything(), ~ sum(is.na(.)))) %>%
  pivot_longer(cols = everything(), names_to = "sample_id", values_to = "NAs")

metadata <- colData(experiment) %>%
  as.data.frame(.) %>%
  rownames_to_column(., var = "sample_id") %>%
  as_tibble(.)

data <- inner_join(metadata, summarized_data, by = 'sample_id') %>%
  group_by(plate) %>%
  mutate(
    outlier = get_zscore_outliers(NAs, sensitivity = 2, sided = "right")
  ) %>%
  ungroup() %>%
  select(
    .,
    -one_of(c("condition", "group", "plate"))
  )

## ---- remove-outliers-and-generate-summarized-experiment ----

f_missingness_cutoff <- z_score_based_missingness_cutoffs['Z > 2']
s_to_exclude <- data[data$outlier == TRUE, ] %>% pull(sample_id)

experiment <- prepare_proteomics_data(
  proteomics_data = global_dat,
  proteomics_meta = meta_data_global,
  type = 'global',
  sample_col_prefix = "G_",
  group_contains_spheroid = TRUE,
  feature_missingness_cutoff = f_missingness_cutoff,
  samples_to_exclude = s_to_exclude,
  imputation_method = "SampMin",
  median_center = TRUE
)
