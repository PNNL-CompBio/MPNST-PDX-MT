##run the DMEA for individual experiments and pooling across samples

##include names of drugs and stuff here
source("plottingParameters.R")

# load packages
if (!require(panSEA)) {
  if (!require("devtools", quietly = TRUE))
    install.packages("devtools")

  devtools::install_github("BelindaBGarana/panSEA", "add-tie-handling")
}

library(plyr);library(dplyr);library(ggplot2);library(synapser);library(panSEA)
library(data.table); library(RColorBrewer)
#source("https://raw.githubusercontent.com/PNNL-CompBio/MPNST_Chr8/refs/heads/main/panSEA_helper_20240913.R")

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

## load inputs
synapser::synLogin()
# RNA-seq data with HGNC symbols (update once batch-corrected)
degs <- read.csv(synapser::synGet('syn69910880')$path)
degs <- degs[degs$hgnc_symbol != "" & !is.na(degs$hgnc_symbol),]
degs <- degs[degs$Drug != "Irinotecan",] # since irinotecan isn't metabolically active
colnames(degs)[ncol(degs)] <- "Gene"
colnames(degs)[4] <- "Log2FC"
degs$padj <- as.numeric(degs$padj)
degs <- degs[degs$padj <= 0.05,]
length(na.omit(unique(degs$Gene)))

# CCLE RNA-seq with HGNC symbols
ccle <- get_CCLE_RNA(hgnc = TRUE)

# reduce to overlapping genes
degs <- degs[degs$Gene %in% colnames(ccle)[2:ncol(ccle)],]
length(na.omit(unique(degs$Gene)))
ccle <- ccle[,c("CCLE_ID", degs$Gene)]
degs$Log2FC <- as.numeric(degs$Log2FC)

# PRISM drug sensitivity (i.e., AUC)
prism <- read.csv(file = paste0("https://raw.github.com/BelindaBGarana/",
                                "DMEA/shiny-app/Inputs/PRISM_drug_mean_AUC_6-23-21.csv"))
prism$X <- NULL

# PRISM drug mechanism of action annotations
gmt.drug <- GSA::GSA.read.gmt(file = paste0(
  "https://raw.github.com/BelindaBGarana/DMEA/shiny-app/",
  "Inputs/MOA_gmt_file_n6_no_special_chars.gmt"))

#########
###run DMEA on individual samples

# compare treated vs. untreated
# DMEA NES < 0: drug mechanisms which may be more toxic to treated samples
synapse_id <- "syn69910841"
timepoints <- na.omit(unique(degs$Timepoint))
mpnst <- na.omit(unique(degs$individualID))
base.path <- getwd()
dmea.path <- file.path(base.path, "DMEA_individual")
dir.create("DMEA_individual")
setwd("DMEA_individual")

dmea.sens <- list()
drug.corr.sens <- list()
for (i in timepoints) {
  cat(i, "\n")

  # filter for time
  time.deg <- degs[degs$Timepoint == i,]

  temp.dmea.sens <- list()
  temp.drug.corr.sens <- list()
  for (m in mpnst) {
    cat(m, "\n")

    # filter for mpnst
    mpnst.deg <- time.deg[time.deg$individualID == m,]

    # prep inputs for each drug
    drug.degs <- list()
    expr <- list()
    all.drugs <- na.omit(unique(mpnst.deg$Drug))
    for (j in all.drugs) {
      cat(j,"(",which(all.drugs == j),"/",length(all.drugs),"drugs)\n")
      temp.deg <- na.omit(dplyr::distinct(
        mpnst.deg[mpnst.deg$Drug == j,
                  c("Gene","Log2FC")]))
      if (any(duplicated(temp.deg$Gene))) {
        cat("resolving duplicates\n")
        temp.deg <- plyr::ddply(temp.deg, .(Gene), summarize,
                                Log2FC = mean(Log2FC))
      }
      temp.deg <- temp.deg[temp.deg$Gene %in% colnames(ccle)[2:ncol(ccle)],]
      if (nrow(temp.deg) > 0) {
        drug.degs[[j]] <- temp.deg
        expr[[j]] <- ccle[,c("CCLE_ID", temp.deg$Gene)]
      }
    }

    # also run analysis for each drug
    drugResults <- panSEA::mDMEA(drug.sensitivity = prism, gmt = gmt.drug,
                                 expression = expr, weights = drug.degs,
                                 types = names(drug.degs), scatter.plots=FALSE)

    temp.dmea.sens[[m]] <- drugResults$compiled.results$result

    # also compile drug correlations
    drug.corr.list <- list()
    for (j in all.drugs) {
      drug.corr.list[[j]] <- drugResults$all.results[[j]]$corr.result
    }
    drug.corr.df <- data.table::rbindlist(drug.corr.list, use.names = TRUE, idcol = "DrugTreatment")
    temp.drug.corr.sens[[m]] <- drug.corr.df
    saveRDS(temp.dmea.sens, file.path(dmea.path,"tempDMEAsens8h.rds"))
    saveRDS(temp.drug.corr.sens, file.path(dmea.path,"tempDMEAcorrSens8h.rds"))
  }
  temp.drug.corr.sens.df <- data.table::rbindlist(temp.drug.corr.sens, use.names=TRUE, idcol="MPNST")
  temp.dmea.sens.df <- data.table::rbindlist(temp.dmea.sens, use.names=TRUE, idcol="MPNST")
  colnames(temp.dmea.sens.df)[2] <- "DrugTreatment"
  temp.dmea.sens.df[temp.dmea.sens.df$DrugTreatment == "Ribociblib",]$DrugTreatment <- "Ribociclib"
  temp.drug.corr.sens.df[temp.drug.corr.sens.df$DrugTreatment == "Ribociblib",]$DrugTreatment <- "Ribociclib"
  drug.corr.sens[[as.character(i)]] <- temp.drug.corr.sens.df
  dmea.sens[[as.character(i)]] <- temp.dmea.sens.df
}
setwd(dmea.path)
drug.corr.sens.df <- data.table::rbindlist(drug.corr.sens, use.names=TRUE, idcol="Time")
dmea.sens.df <- data.table::rbindlist(dmea.sens, use.names=TRUE, idcol="Time")
write.csv(dmea.sens.df,"DMEA_results.csv", row.names = FALSE)
write.csv(drug.corr.sens.df,"DMEA_correlations.csv", row.names = FALSE)
synapser::synStore(synapser::File("DMEA_results.csv", parent=synapse_id))
synapser::synStore(synapser::File("DMEA_correlations.csv", parent=synapse_id))


#################################
## run DMEA on pooled samples

# synapser::synLogin()
 # RNA-seq data with HGNC symbols
 degs <- read.csv(synapser::synGet('syn64397743')$path) # 1312039 genes
 degs <- degs[degs$hgnc_symbol != "" & !is.na(degs$hgnc_symbol),]
 degs <- degs[degs$Drug != "Irinotecan",] # since irinotecan isn't metabolically active
 colnames(degs)[ncol(degs)] <- "Gene"
 colnames(degs)[4] <- "Log2FC"
 degs$padj <- as.numeric(degs$padj)
 degs <- degs[degs$padj <= 0.05,]
 length(na.omit(unique(degs$Gene))) # 11194

drug.info <- list("CDK inhibitor" = c("Palbociclib","Ribociclib"),
                  "Chemotherapy" = "Trabectedin",
                  "DNA alkylating agent" = "Ifosfamide",
                  "DNA methyltransferase inhibitor" = "Decitabine",
                  "HDAC inhibitor" = "Vorinostat",
                  "MEK inhibitor" = c("Mirdametinib","Selumetinib","Trametinib"),
                  "MET inhibitor" = "Capmatinib",
                  "PARP inhibitor" = "Olaparib",
                  "SHP2 inhibitor" = c("RMC4630","TNO155"),
                  "topoisomerase inhibitor" = c("Doxorubicin", "Irinotecan"),
                  "TEADi" = "Verteporfin")
# compare treated vs. untreated
# DMEA NES < 0: drug mechanisms which may be more toxic to treated samples
synapse_id <- "syn64618500"
timepoints <- na.omit(unique(degs$Timepoint))
mpnst <- na.omit(unique(degs$individualID))
base.path <- getwd()
dmea.path <- file.path(base.path, "DMEA")
dir.create("DMEA")
setwd("DMEA")

dmea.sens <- list()
drug.corr.sens <- list()
for (i in timepoints) {
  cat(i, "\n")

  # filter for time
  time.deg <- degs[degs$Timepoint == i,]

  temp.dmea.sens <- list()
  temp.drug.corr.sens <- list()

  # prep inputs for each drug
  drug.degs <- list()
  expr <- list()
  all.drugs <- na.omit(unique(time.deg$Drug))
  for (j in all.drugs) {
    cat(j,"(",which(all.drugs == j),"/",length(all.drugs),"drugs)\n")
    temp.deg <- na.omit(dplyr::distinct(
      time.deg[time.deg$Drug == j,
               c("Gene","Log2FC")]))
    if (any(duplicated(temp.deg$Gene))) {
      cat("resolving duplicates\n")
      temp.deg <- plyr::ddply(temp.deg, .(Gene), summarize,
                              Log2FC = mean(Log2FC))
    }
    temp.deg <- temp.deg[temp.deg$Gene %in% colnames(ccle)[2:ncol(ccle)],]
    if (nrow(temp.deg) > 0) {
      drug.degs[[j]] <- temp.deg
      expr[[j]] <- ccle[,c("CCLE_ID", temp.deg$Gene)]
    }
  }

  # also run analysis for each drug
  drugResults <- panSEA::mDMEA(drug.sensitivity = prism, gmt = gmt.drug,
                               expression = expr, weights = drug.degs,
                               types = names(drug.degs), scatter.plots=FALSE)

  temp.dmea.sens.df <- drugResults$compiled.results$result

  # also compile drug correlations
  drug.corr.list <- list()
  for (j in all.drugs) {
    drug.corr.list[[j]] <- drugResults$all.results[[j]]$corr.result
  }
  temp.drug.corr.sens.df <- data.table::rbindlist(drug.corr.list, use.names = TRUE, idcol = "DrugTreatment")
  colnames(temp.dmea.sens.df)[1] <- "DrugTreatment"
  drug.corr.sens[[as.character(i)]] <- temp.drug.corr.sens.df
  dmea.sens[[as.character(i)]] <- temp.dmea.sens.df
}
setwd(dmea.path)
drug.corr.sens.df <- data.table::rbindlist(drug.corr.sens, use.names=TRUE, idcol="Time")
dmea.sens.df <- data.table::rbindlist(dmea.sens, use.names=TRUE, idcol="Time")
write.csv(dmea.sens.df,"DMEA_results.csv", row.names = FALSE)
write.csv(drug.corr.sens.df,"DMEA_correlations.csv", row.names = FALSE)
synapser::synStore(synapser::File("DMEA_results.csv", parent=synapse_id))
synapser::synStore(synapser::File("DMEA_correlations.csv", parent=synapse_id))
