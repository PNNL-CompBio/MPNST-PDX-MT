##run GSEA HALLMARK analysis on bulk and individual samples

# load packages
if (!require(panSEA)) {
  if (!require("devtools", quietly = TRUE))
    install.packages("devtools")

  devtools::install_github("PNNL-CompBio/panSEA", "add-tie-handling")
}

library(plyr);library(dplyr);library(ggplot2); library(patchwork)
library(synapser);library(panSEA); library(data.table);
library(ComplexHeatmap)

# RNA-seq data with HGNC symbols
diffexp <- read.csv(synapser::synGet('syn64397743')$path)
diffexp <- na.omit(diffexp[diffexp$hgnc_symbol!="",
                           c("hgnc_symbol", "log2FoldChange", "Drug",
                             "Timepoint")])
diffexp$log2FoldChange <- as.numeric(diffexp$log2FoldChange)
diffexp <- diffexp[diffexp$Drug != "Irinotecan",] # since irinotecan isn't metabolically active

# get gene sets: Hallmark (H)
gmt.Hallmark <- as.data.frame(msigdbr::msigdbr(collection = "H"))
gmt.Hallmark <- DMEA::as_gmt(gmt.Hallmark, element.names="gene_symbol", set.names="gs_name") # 196 gene sets
saveRDS(gmt.Hallmark,"gmt_Hallmark.rds")

synapse_id <- "syn65916546"
timepoints <- na.omit(unique(diffexp$Timepoint))
base.path <- getwd()
gsea.path <- file.path(base.path, "GSEA_Hallmark")
dir.create("GSEA_Hallmark")

gsea <- list()
for (i in timepoints) {
  # filter for time
  time.deg <- diffexp[diffexp$Timepoint == i,]

  # prep inputs for each drug
  drug.diffexp <- list()
  hall.gmts <- list()
  all.drugs <- na.omit(unique(time.deg$Drug))
  for (j in all.drugs) {
    temp.deg <- na.omit(dplyr::distinct(
      time.deg[time.deg$Drug == j,c("hgnc_symbol","log2FoldChange")]))
    colnames(temp.deg) <- c("Gene", "Log2FC")
    if (any(duplicated(temp.deg$Gene))) {
      temp.deg <- plyr::ddply(temp.deg, .(Gene), summarize,
                              Log2FC = mean(Log2FC))
    }
    drug.diffexp[[j]] <- temp.deg
    hall.gmts[[j]] <- gmt.Hallmark
  }

  # also run analysis for each drug
  hallResults <- panSEA::mGSEA(drug.diffexp, gmt=hall.gmts,
                               types = names(drug.diffexp))
  temp.gsea.df <- hallResults$compiled.results$result
  temp.gsea.df$Collection <- "Hallmark"
  colnames(temp.gsea.df)[1] <- "DrugTreatment"
  gsea[[as.character(i)]] <- temp.gsea.df
}
setwd(gsea.path)
gsea.df <- data.table::rbindlist(gsea, use.names=TRUE, idcol="Time")
write.csv(gsea.df,"GSEA_results_Hallmark.csv", row.names = FALSE)
synapser::synStore(synapser::File("GSEA_results_Hallmark.csv", parent=synapse_id))


#######now run individual results

# RNA-seq data with HGNC symbols
diffexp <- read.csv(synapser::synGet('syn69910880')$path)
diffexp <- na.omit(diffexp[diffexp$hgnc_symbol!="",
                           c("hgnc_symbol", "log2FoldChange", "Drug",
                             "Timepoint","individualID")])
diffexp$log2FoldChange <- as.numeric(diffexp$log2FoldChange)
diffexp <- diffexp[diffexp$Drug != "Irinotecan",] # since irinotecan isn't metabolically active
# GSEA NES > 0: pathways up-regulated in treated vs. untreated samples
synapse_id <- "syn69910857"
timepoints <- na.omit(unique(diffexp$Timepoint))
base.path <- getwd()
gsea.path <- file.path(base.path, "GSEA_Hallmark_individual")
dir.create("GSEA_Hallmark_individual")
mpnst<-unique(diffexp$individualID) ##ADDed by SG
gsea <- list()
for (i in timepoints) {
  # filter for time
  time.deg <- diffexp[diffexp$Timepoint == i,]
  temp.gsea <- list()
  for (m in mpnst) {
    # filter for mpnst
    mpnst.deg <- time.deg[time.deg$individualID == m,]

    # prep inputs for each drug
    drug.diffexp <- list()
    kegg.gmts <- list()
    hall.gmts <- list()
    all.drugs <- na.omit(unique(mpnst.deg$Drug))
    for (j in all.drugs) {
      temp.deg <- na.omit(dplyr::distinct(
        mpnst.deg[mpnst.deg$Drug == j,c("hgnc_symbol","log2FoldChange")]))
      colnames(temp.deg) <- c("Gene", "Log2FC")
      if (any(duplicated(temp.deg$Gene))) {
        temp.deg <- plyr::ddply(temp.deg, .(Gene), summarize,
                                Log2FC = mean(Log2FC))
      }
      drug.diffexp[[j]] <- temp.deg
      kegg.gmts[[j]] <- gmt.Hallmark
    }

    # also run analysis for each drug
    keggResults <- panSEA::mGSEA(drug.diffexp, gmt=kegg.gmts,
                                 types = names(drug.diffexp))
    temp.gsea.kegg <- keggResults$compiled.results$result
    temp.gsea.kegg$Collection <- "Hallmark"
    temp.gsea[[m]] <- temp.gsea.kegg
  }
  temp.gsea.df <- data.table::rbindlist(temp.gsea, use.names=TRUE, idcol="MPNST")
  colnames(temp.gsea.df)[2] <- "DrugTreatment"
  gsea[[as.character(i)]] <- temp.gsea.df
}
setwd(gsea.path)
gsea.df <- data.table::rbindlist(gsea, use.names=TRUE, idcol="Time")
write.csv(gsea.df,"GSEA_results_Hallmark.csv", row.names = FALSE)
synapser::synStore(synapser::File("GSEA_results_Hallmark.csv", parent=synapse_id))

