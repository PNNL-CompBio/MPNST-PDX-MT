# MEK+HDACi
library(plyr);library(dplyr);library(synapser);
library(ggplot2);library(ComplexHeatmap);library(RColorBrewer); library(ggvenn)
#source("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/helperScripts/compile_mDEG.R")
library(panSEA)

setwd("~/OneDrive - PNNL/Documents/MPNST")
dir.create("MEKi")
setwd("MEKi")
dir.create("Treated")
setwd("Treated")
dir.create("GSE183307_and_GSE262030")
setwd("GSE183307_and_GSE262030")
dir.create("24h")
setwd("24h")
dir.create("padj0.05")
setwd("padj0.05")
dir.create("min2PerGroup")
setwd("min2PerGroup")

#### mDMEA with sig filter ####
exps <- c("GSE183307", "GSE262030")
weights <- list()
for (i in exps) {
  temp.sig <- read.csv(paste0(i, "_diffexp.csv"))
  weights[[i]] <- na.omit(temp.sig[temp.sig$adj.P.Val <= 0.05,])
}
#source("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/helperScripts/compile_mDEG.R")
#degResults <- compile_mDEG(weights)

venn.list <- list()
for (i in exps) {
  venn.list[[i]] <- na.omit(weights[[i]])$Gene
}
ggvenn::ggvenn(venn.list, show_percentage = FALSE, text_size=5, set_name_size=5)
ggplot2::ggsave("diffexp_vennDiagram.pdf", width=5, height=5)

# merge for correlation and scatter plot
deg <- data.frame()
for (i in exps) {
  deg1 <- weights[[i]]
  deg1$Source <- i
  deg <- rbind(deg, deg1)
}
combined.degs <- na.omit(reshape2::dcast(deg, Gene ~ Source, value.var="Log2FC")) # 578
pearson.test <- cor.test(combined.degs[,2], combined.degs[,3])
pearson.p <- pearson.test$p.value
pearson.est <- pearson.test$estimate
stats_pearson <- substitute(
  r == est * "," ~ ~"p" ~ "=" ~ p,
  list(
    est = as.numeric(format(pearson.est, digits = 3)),
    p = format(pearson.p, digits = 3)
  )
)
ggplot2::ggplot(combined.degs, aes(x=GSE183307, y=GSE262030)) + geom_point() + 
  theme_minimal(base_size=12) + ggplot2::geom_smooth(method = "lm", size = 1.5, 
                       linetype = "solid", color = "blue", se = TRUE, na.rm = TRUE) +
  ggplot2::geom_text(x = max(combined.degs$GSE183307), y = min(combined.degs$GSE262030),
                     vjust = "inward", hjust = "inward", colour = "blue", parse = TRUE,
    label = as.character(as.expression(stats_pearson)), size = 8) + 
  ggtitle(paste0("Log2 Fold-change (N = ",nrow(combined.degs),")")) +
  geom_hline(yintercept=0) + geom_vline(xintercept=0)
ggplot2::ggsave("diffexp_scatterPlot_Pearson.pdf", width = 5, height=5)

spearman.test <- cor.test(combined.degs[,2], combined.degs[,3], method="spearman")
spearman.p <- spearman.test$p.value
spearman.est <- spearman.test$estimate
stats_spearman <- substitute(
  rho == est * "," ~ ~"p" ~ "=" ~ p,
  list(
    est = as.numeric(format(spearman.est, digits = 3)),
    p = format(spearman.p, digits = 3)
  )
)
ggplot2::ggplot(combined.degs, aes(x=rank(GSE183307), y=rank(GSE262030))) + geom_point() + 
  theme_classic(base_size=12) + ggplot2::geom_smooth(method = "lm", size = 1.5, 
                                                     linetype = "solid", color = "blue", se = TRUE, na.rm = TRUE) +
  ggplot2::geom_text(x = max(rank(combined.degs$GSE183307)), y = min(rank(combined.degs$GSE262030)),
                     vjust = "inward", hjust = "inward", colour = "blue", parse = TRUE,
                     label = as.character(as.expression(stats_spearman)), size = 8) + 
  ggtitle(paste0("Log2 Fold-change (N = ",nrow(combined.degs),")"))
ggplot2::ggsave("diffexp_scatterPlot_Spearman.pdf", width = 5, height=5)
nrow(combined.degs[combined.degs$GSE183307/combined.degs$GSE262030 > 0,]) # 1034 (94%) are in same direction
nrow(combined.degs[combined.degs$GSE183307/combined.degs$GSE262030 < 0,]) # 66 (6%) are in opposite directions

dmeaResults <- panSEA::mDMEA(weights=weights, types=exps)

source("https://raw.githubusercontent.com/PNNL-CompBio/MPNST_Chr8/refs/heads/main/panSEA_helper_20240913.R")
temp.files <- list()
for (i in exps) {
  temp.files[[i]] <- extract_DMEA_files(dmeaResults, i)
}
temp.files[["Compiled"]] <- list("DMEA_results.csv" = dmeaResults$compiled.results$results,
                                 "Compiled_DMEA_results.csv" = dmeaResults$compiled.results$mean.results,
                                 "DMEA_venn_diagram.pdf" = dmeaResults$compiled.results$venn.diagram,
                                 "DMEA_dot_plot.pdf" = dmeaResults$compiled.results$dot.plot,
                                 "DMEA_correlations.csv" = dmeaResults$compiled.results$corr,
                                 "DMEA_correlation_matrix.pdf" =dmeaResults$compiled.results$corr.matrix)
save_to_synapse_v2(temp.files)

# redo dot plot
# reduce plot data down to top results
sig.DMEA.df <- dmeaResults$compiled.results$mean.results[dmeaResults$compiled.results$mean.results$N_sig > 0, ]
n.dot.sets <- 10
if (nrow(sig.DMEA.df) >= n.dot.sets) {
  top.DMEA.df <- 
    sig.DMEA.df %>% dplyr::slice_max(abs(mean_NES), n = n.dot.sets)
} else {
  top.DMEA.df <- 
    dmeaResults$compiled.results$mean.results %>% dplyr::slice_max(abs(mean_NES), n = n.dot.sets)
}
dot.df <- dmeaResults$compiled.results$results[dmeaResults$compiled.results$results$Drug_set %in% top.DMEA.df$Drug_set, ]
dot.df$minusLogFDR <- -log10(1E-4)
dot.df[dot.df$FDR_q_value!=0,]$minusLogFDR <- -log10(dot.df[dot.df$FDR_q_value!=0,]$FDR_q_value)
mean.DMEA.df <- dplyr::arrange(dmeaResults$compiled.results$mean.results, desc(mean_NES))

library(ggplot2)
maxAbsNES <- max(abs(dot.df$NES))
dot.plot <- ggplot2::ggplot(
  dot.df,
  ggplot2::aes(
    x = type, y = Drug_set, color = NES,
    size = minusLogFDR
  )
) +
  ggplot2::geom_point() +
  ggplot2::scale_y_discrete(limits = mean.DMEA.df[
    mean.DMEA.df$Drug_set %in% top.DMEA.df$Drug_set, ]$Drug_set) +
  scale_color_gradient2(low="blue",high="red", mid="grey", 
                        limits=c(-maxAbsNES, maxAbsNES)) +
  theme_classic(base_size=12) + 
  theme(axis.title=element_blank(), 
        axis.text.x = element_text(angle=45, vjust=1, hjust=1)) +
  ggplot2::labs(color = "NES", size = "-log(FDR)") + 
  geom_point(data = subset(dot.df, sig), col = "black", stroke = 1.5, shape = 21)
dot.plot
ggplot2::ggsave("Compiled_DMEA_dotPlot.pdf", width=4, height=3.75)
ggplot2::ggsave("Compiled_DMEA_dotPlot_taller.pdf", width=4, height=4)

# generate list of significant results
sigList <- list()
for (i in exps) {
  sigList[[i]] <- na.omit(unique(dmeaResults$compiled.results$results[dmeaResults$compiled.results$results$type == i & dmeaResults$compiled.results$results$sig,]$Drug_set))
}

# convert to matrix
ggvenn::ggvenn(sigList, show_percentage = FALSE, set_name_size = 5, text_size = 5)
ggplot2::ggsave("MOA_sensitivity_vennDiagram.pdf", width=5, height=5)

# merge for correlation and scatter plot
dmeaDF <- data.frame()
for (i in exps) {
  dmeaDF1 <- na.omit(unique(dmeaResults$compiled.results$results[dmeaResults$compiled.results$results$type == i & dmeaResults$compiled.results$results$sig,]))
  dmeaDF1$Source <- i
  dmeaDF <- rbind(dmeaDF, dmeaDF1)
}
combined.degs <- na.omit(reshape2::dcast(dmeaDF, Drug_set ~ type, value.var="NES")) # 578
pearson.test <- cor.test(combined.degs[,2], combined.degs[,3])
pearson.p <- pearson.test$p.value
pearson.est <- pearson.test$estimate
stats_pearson <- substitute(
  r == est * "," ~ ~"p" ~ "=" ~ p,
  list(
    est = as.numeric(format(pearson.est, digits = 3)),
    p = format(pearson.p, digits = 3)
  )
)
ggplot2::ggplot(combined.degs, aes(x=GSE183307, y=GSE262030)) + geom_point() + 
  theme_minimal(base_size=12) + ggplot2::geom_smooth(method = "lm", size = 1.5, 
                                                     linetype = "solid", color = "blue", se = TRUE, na.rm = TRUE) +
  ggplot2::geom_text(x = max(combined.degs$GSE183307), y = min(combined.degs$GSE262030),
                     vjust = "inward", hjust = "inward", colour = "blue", parse = TRUE,
                     label = as.character(as.expression(stats_pearson)), size = 8) + 
  ggtitle(paste0("NES (N = ",nrow(combined.degs),")")) +
  geom_hline(yintercept=0) + geom_vline(xintercept=0)
ggplot2::ggsave("DMEA_sensitivity_scatterPlot_Pearson.pdf", width = 5, height=5)
ggplot2::ggplot(combined.degs, aes(x=GSE183307, y=GSE262030)) + geom_point() + 
  theme_classic(base_size=12) + ggplot2::geom_smooth(method = "lm", size = 1.5, 
                                                     linetype = "solid", color = "blue", se = TRUE, na.rm = TRUE) +
  ggplot2::geom_text(x = max(combined.degs$GSE183307), y = min(combined.degs$GSE262030),
                     vjust = "inward", hjust = "inward", colour = "blue", parse = TRUE,
                     label = as.character(as.expression(stats_pearson)), size = 8) + 
  ggtitle(paste0("NES (N = ",nrow(combined.degs),")"))
ggplot2::ggsave("DMEA_sensitivity_scatterPlot_Pearson_v2.pdf", width = 5, height=5)

spearman.test <- cor.test(combined.degs[,2], combined.degs[,3], method="spearman")
spearman.p <- spearman.test$p.value
spearman.est <- spearman.test$estimate
stats_spearman <- substitute(
  rho == est * "," ~ ~"p" ~ "=" ~ p,
  list(
    est = as.numeric(format(spearman.est, digits = 3)),
    p = format(spearman.p, digits = 3)
  )
)
ggplot2::ggplot(combined.degs, aes(x=rank(GSE183307), y=rank(GSE262030))) + geom_point() + 
  theme_classic(base_size=12) + ggplot2::geom_smooth(method = "lm", size = 1.5, 
                                                     linetype = "solid", color = "blue", se = TRUE, na.rm = TRUE) +
  ggplot2::geom_text(x = max(rank(combined.degs$GSE183307)), y = min(rank(combined.degs$GSE262030)),
                     vjust = "inward", hjust = "inward", colour = "blue", parse = TRUE,
                     label = as.character(as.expression(stats_spearman)), size = 8) + 
  ggtitle(paste0("NES (N = ",nrow(combined.degs),")"))
ggplot2::ggsave("DMEA_sensitivity_scatterPlot_Spearman.pdf", width = 5, height=5)
nrow(combined.degs[combined.degs$GSE183307/combined.degs$GSE262030 > 0,]) # 6 (100%) are in same direction
nrow(combined.degs[combined.degs$GSE183307/combined.degs$GSE262030 < 0,]) # 0 (0%) are in opposite directions

#### Drug correlation bar plot with sig filter ####
sig.moas <- unique(dmeaResults$compiled.results$mean.results[dmeaResults$compiled.results$mean.results$N_sig>0,]$Drug_set)
drug.corr <- list()
for (i in exps) {
  drug.corr[[i]] <- dmeaResults$all.results[[i]]$corr.result
}
drug.corr <- data.table::rbindlist(drug.corr, use.names=TRUE, idcol="Source")
drug.info <- read.csv("https://raw.githubusercontent.com/BelindaBGarana/DMEA/refs/heads/shiny-app/Inputs/PRISM_secondary-screen-replicate-treatment-info.csv")
maxAbsEst <- max(na.omit(abs(drug.corr$Pearson.est)))
#maxAbsEst <- 1
drug.corr <- na.omit(drug.corr)
if (any(drug.corr$Pearson.q == 0)) {
  drug.corr[drug.corr$Pearson.q == 0,]$Pearson.q <- 0.0001
}
drug.corr$minusLogP <- -log(drug.corr[,"Pearson.p"], base = 10)
drug.corr$minusLogFDR <- -log(drug.corr[,"Pearson.q"], base = 10)
drug.info$Drug <- gsub("-",".",drug.info$name)
drug.info$Drug <- gsub("[(]",".",drug.info$Drug)
drug.info$Drug <- gsub("[)]",".",drug.info$Drug)
drug.info$Drug <- gsub("[+]",".",drug.info$Drug)
red.drug.info <- dplyr::distinct(drug.info[,c("Drug","name","moa","target")])
colnames(drug.corr)[2] <- "Drug"
drug.corr.wInfo <- merge(red.drug.info, drug.corr, by="Drug", all.y = TRUE)
drug.corr.wInfo$Mechanism <- "Other"
drug.corr.wInfo[drug.corr.wInfo$moa %in% sig.moas,]$Mechanism <- drug.corr.wInfo[drug.corr.wInfo$moa %in% sig.moas,]$moa

maxAbsEst <- max(na.omit(abs(drug.corr.wInfo$Pearson.est)))
maxLogFDR <- ceiling(max(na.omit(drug.corr.wInfo$minusLogFDR)))
library(patchwork); library(ggplot2)
moaOrder <- dmeaResults$compiled.results$mean.results[order(dmeaResults$compiled.results$mean.results$mean_NES),]$Drug_set # 24
#drug.corr.wInfo$Mechanism <- factor(drug.corr.wInfo$Mechanism, levels = c(moaOrder, "Other"))
drug.corr.wInfo$Significant <- FALSE
drug.corr.wInfo[drug.corr.wInfo$Pearson.q <= 0.05,]$Significant <- TRUE
drug.corr.wInfo <- drug.corr.wInfo[!is.na(drug.corr.wInfo$name),]

MOAsInTop50 <- c()
for (i in exps) {
  temp.dot.df <- drug.corr.wInfo[drug.corr.wInfo$Significant & drug.corr.wInfo$Source == i,]
  temp.dot.df <- temp.dot.df %>% slice_max(abs(Pearson.est),n=50)
  MOAsInTop50 <- unique(c(MOAsInTop50, temp.dot.df$Mechanism))
}

# take care of cases where there are multiple mechanisms
for (i in 1:length(moaOrder)) {
  if (nrow(drug.corr.wInfo[drug.corr.wInfo$Mechanism == "Other" & 
                           grepl(moaOrder[i],drug.corr.wInfo$moa),]) > 0) {
    drug.corr.wInfo[drug.corr.wInfo$Mechanism == "Other" & 
                      grepl(moaOrder[i],drug.corr.wInfo$moa),]$Mechanism <- moaOrder[i]
  }
}

# repeat in case of any change
MOAsInTop50 <- c() # comes out to 17
for (i in exps) {
  temp.dot.df <- drug.corr.wInfo[drug.corr.wInfo$Significant & drug.corr.wInfo$Source == i,]
  temp.dot.df <- temp.dot.df %>% slice_max(abs(Pearson.est),n=50)
  MOAsInTop50 <- unique(c(MOAsInTop50, temp.dot.df$Mechanism))
}
moaOrder3 <- moaOrder[moaOrder %in% MOAsInTop50] # 28

maxAbsEst <- max(na.omit(abs(drug.corr.wInfo$Pearson.est)))
maxLogFDR <- ceiling(max(na.omit(drug.corr.wInfo$minusLogFDR)))
library(patchwork); library(ggplot2)
drug.corr.wInfo$Significant <- FALSE
drug.corr.wInfo[drug.corr.wInfo$Pearson.q <= 0.05,]$Significant <- TRUE

library(plyr);library(dplyr)

gsea.dot.plots <- NULL
for (i in exps) {
  nGenes <- length(unique(na.omit(drug.corr.wInfo[drug.corr.wInfo$Source == i,])$Drug))
  nCorrGenes <- length(unique(na.omit(drug.corr.wInfo[drug.corr.wInfo$Significant & 
                                                        drug.corr.wInfo$Source == i,])$Drug))
  
  temp.dot.df <- drug.corr.wInfo[drug.corr.wInfo$Significant & drug.corr.wInfo$Source == i,]
  temp.dot.df <- temp.dot.df %>% slice_max(abs(Pearson.est),n=50)
  topGenes <- temp.dot.df$name
  nTopGenes <- length(topGenes)
  geneOrder <- temp.dot.df[order(temp.dot.df$Pearson.est),]$name
  
  plot.annot <- paste0(i, "\n(", nCorrGenes, " / ", nGenes, " drugs correlated)")
  
  dot.plot <- ggplot2::ggplot(temp.dot.df,
                              ggplot2::aes(
                                x = name, y = Pearson.est, fill = Mechanism
                              )
  ) + scale_size(limits=c(0,maxLogFDR), range = c(0.5,4)) +
    ggplot2::geom_col() +
    ggplot2::scale_x_discrete(limits = geneOrder) +
    theme_classic() + scale_fill_manual(breaks=c(moaOrder3,"Other"), 
                                        values = c(grDevices::colorRampPalette(
                                          RColorBrewer::brewer.pal(12, "Set3"))(length(moaOrder3)))) +
    ggplot2::labs(
      #x = i,
      y = "Pearson r",
      fill = "Drug Mechanism"
    ) + theme(axis.title.x = element_blank()) +
    #theme(axis.title.y=element_text(face="bold", size=16), 
    #axis.text.y = element_blank(),
    #axis.ticks.y = element_blank(),
    #axis.title.x=element_blank()) +
    theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1)) +
    theme(legend.direction = "horizontal", legend.position = "bottom")
  
  cat(plot.annot, "\n")
  cat(i, "with", nTopGenes, "top Drugs\n")
  dot.plot <- dot.plot + ggtitle(plot.annot) + theme(plot.title = element_text(hjust = 0.5, face="bold", size=16))
  dot.plot
  if (is.null(gsea.dot.plots)) {
    gsea.dot.plots <- dot.plot
  } else {
    gsea.dot.plots <- gsea.dot.plots / dot.plot
  } 
}
source("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/MPNST/Chr8/MPNST_Chr8_manuscript/Figure_3_Kinase/guides_build_mod.R")
gsea.dot.plots2 <- gsea.dot.plots / plot_spacer() + plot_layout(guides = 'collect')
gsea.dot.plots2
ggplot2::ggsave("CorrelatedDrugs_barPlot_patchworkSource_moaFill_sliceMaxAbsPearson50_oppositeMOAorder.pdf", gsea.dot.plots, width=10, height=8)
