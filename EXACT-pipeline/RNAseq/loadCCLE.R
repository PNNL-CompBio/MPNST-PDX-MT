##load drug data for DMEA
#moving to separate file for use by other tools


##function to lod mRNA

# load helper functions
#source("https://github.com/PNNL-CompBio/coderdata/blob/main/build/lincs/05-LINCS_perturbations.R")
get_CCLE_RNA <- function(hgnc = FALSE) {
  message("Loading adherent CCLE RNA-seq data version 19Q4")
  if (hgnc) {
    if (!file.exists("Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_approved_symbols_only_1-200.Rbin")) {
      download.file(
        paste0(
          "https://raw.github.com/BelindaBGarana/DMEA/shiny-app/Inputs/",
          "Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_approved_symbols_only_1-200.Rbin"
        ),
        destfile =
          "Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_approved_symbols_only_1-200.Rbin"
      )
    }

    if (!file.exists("Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_approved_symbols_only_201-327.Rbin")) {
      download.file(
        paste0(
          "https://raw.github.com/BelindaBGarana/DMEA/shiny-app/Inputs/",
          "Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_approved_symbols_only_201-327.Rbin"
        ),
        destfile =
          "Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_approved_symbols_only_201-327.Rbin"
      )
    }

    RNA.first200 <- readRDS("Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_approved_symbols_only_1-200.Rbin")
    RNA.rest <- readRDS("Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_approved_symbols_only_201-327.Rbin")
    RNA.df <- rbind(RNA.first200, RNA.rest)
  } else {
    if (!file.exists("Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_1-200.Rbin")) {
      download.file(
        paste0(
          "https://raw.github.com/BelindaBGarana/DMEA/shiny-app/Inputs/",
          "Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_1-200.Rbin"
        ),
        destfile =
          "Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_1-200.Rbin"
      )
    }

    if (!file.exists("Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_201-327.Rbin")) {
      download.file(
        paste0(
          "https://raw.github.com/BelindaBGarana/DMEA/shiny-app/Inputs/",
          "Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_201-327.Rbin"
        ),
        destfile =
          "Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_201-327.Rbin"
      )
    }

    load("Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_1-200.Rbin")
    load("Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_201-327.Rbin")
    RNA.df <- rbind(RNA.first200, RNA.rest)
  }
  return(RNA.df)
}



# PRISM drug sensitivity (i.e., AUC)
prism <- read.csv(file = paste0("https://raw.github.com/BelindaBGarana/",
                                "DMEA/shiny-app/Inputs/PRISM_drug_mean_AUC_6-23-21.csv"))
prism$X <- NULL

# PRISM drug mechanism of action annotations
gmt.drug <- GSA::GSA.read.gmt(file = paste0(
  "https://raw.github.com/BelindaBGarana/DMEA/shiny-app/",
  "Inputs/MOA_gmt_file_n6_no_special_chars.gmt"))
