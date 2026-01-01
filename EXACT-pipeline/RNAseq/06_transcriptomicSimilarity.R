# compile transcriptomic results looking for similarity
rm(list=ls())
library(plyr);library(dplyr);library(synapser);library(data.table);library(ggplot2);library(patchwork)
synapser::synLogin()
#setwd("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/MPNST-PDX-MT/RNAseq/synapse_based_code")
dir.create("similarity")
setwd("similarity")


DMEAdotPlots <- function(dot.df,  title="", fname, min.width = 5, min.height = 4,
                         maxAbsNES = max(abs(dot.df[,color])),
                         maxLogFDR = max(dot.df[,size]), x.order=sort(unique(dot.df$DrugTreatment)),
                         colorLab = "NES", width=NULL, height=NULL) {
  dot.plot <- ggplot2::ggplot(dot.df, ggplot2::aes(x = DrugTreatment, y = reorder(Drug_set, NES),
                                                   color = NES, size = minusLogFDR)) +
    facet_grid(Time ~ MPNST) + ggplot2::geom_point() +
    #ggplot2::scale_y_discrete(limits = y.order) +
    ggplot2::scale_x_discrete(limits = x.order) +
    scale_color_gradient2(low="blue",high="red", mid="grey",
                          limits=c(-maxAbsNES, maxAbsNES)) +
    theme_classic(base_size=12) + scale_size_continuous(limits=c(0,maxLogFDR)) +
    ggplot2::labs(color = "NES", size = "-log(FDR)", title=title) +
    theme(axis.title=element_blank(), plot.title=element_text(hjust=0.5, face="bold"),
          axis.text.x = element_text(angle=45, vjust=1, hjust=1)) +
    geom_point(data = subset(dot.df, sig), col = "black", stroke = 1.5, shape = 21)
  customHeight <- ifelse(is.null(height), ifelse(length(unique(dot.df$Drug_set))/5 < min.height, min.height, length(unique(dot.df$Drug_set))/5), height)
  customWidth <- ifelse(is.null(width), min.width, width)
  ggplot2::ggsave(paste0(fname,".pdf"), dot.plot, width=customWidth, height=customHeight)
  ggplot2::ggsave(paste0(fname,"_taller.pdf"), dot.plot, width=customWidth, height=(customHeight+1))
  ggplot2::ggsave(paste0(fname,"_widerTaller.pdf"), dot.plot, width=(customWidth+2), height=(customHeight+2))
  ggplot2::ggsave(paste0(fname,"_evenWiderTaller.pdf"), dot.plot, width=(customWidth+3), height=(customHeight+3))
  ggplot2::ggsave(paste0(fname,"_evenTaller.pdf"), dot.plot, width=customWidth, height=(customHeight+3))
  ggplot2::ggsave(paste0(fname,"_evenWider.pdf"), dot.plot, width=(customWidth+3), height=(customHeight+2))

  if (all(grepl(" inhibitor", unique(dot.df$Drug_set)))) {
    dot.df$Inhibitor <- sub(" inhibitor", "", dot.df$Drug_set)
    dot.plot <- ggplot2::ggplot(
      dot.df,
      ggplot2::aes(
        x = DrugTreatment, y = reorder(Inhibitor, NES), color = NES,
        size = minusLogFDR
      )
    ) + facet_grid(Time ~ MPNST) +
      ggplot2::geom_point() +
      #ggplot2::scale_y_discrete(limits = sigOrder2) +
      ggplot2::scale_x_discrete(limits = x.order) +
      scale_color_gradient2(low="blue",high="red", mid="grey",
                            limits=c(-maxAbsNES, maxAbsNES)) +
      theme_classic(base_size=12) + scale_size_continuous(limits=c(0,maxLogFDR)) +
      ggplot2::labs(color = "NES", size = "-log(FDR)",title=title) +
      theme(axis.title=element_blank(), plot.title=element_text(hjust=0.5, face="bold"),
            axis.text.x = element_text(angle=45, vjust=1, hjust=1)) +
      geom_point(data = subset(dot.df, sig), col = "black", stroke = 1.5, shape = 21)
    ggplot2::ggsave(paste0(fname,"_inhibitors.pdf"), dot.plot, width=(min.width-1), height=customHeight)
    ggplot2::ggsave(paste0(fname,"_inhibitors_wider.pdf"), dot.plot, width=min.width, height=customHeight)
    ggplot2::ggsave(paste0(fname,"_inhibitors_evenWider.pdf"), dot.plot, width=(min.width+1), height=customHeight)
  }
  # else if (customHeight<49) {
  #   ggplot2::ggsave(paste0(fname,"_wider.pdf"), dot.plot, width=(min.width+1), height=customHeight)
  #   ggplot2::ggsave(paste0(fname,"_widerTaller.pdf"), dot.plot, width=(min.width+1), height=(customHeight+1))
  # }
}


# load results
diffexp <- read.csv(synapser::synGet('syn64397743')$path)
diffexp$log2FoldChange <- as.numeric(diffexp$log2FoldChange)
diffexp$padj <- as.numeric(diffexp$padj)
diffexp$Timepoint <- diffexp$Time
diffexp$Timepoint <- factor(diffexp$Timepoint, levels=c("8h", "24h"))
diffexp$feature <- diffexp$hgnc_symbol
diffexp$minusLogFDR <- -log10(diffexp$padj)
diffexp <- diffexp[!is.na(diffexp$hgnc_symbol) & diffexp$hgnc_symbol != "",]
diffexp[diffexp$Drug=="Trabectidin",]$Drug <- "Trabectedin"

sens <- read.csv(synapser::synGet("syn65672258")$path)
sens$sig <- FALSE
sens[sens$p_value <= 0.05 & sens$FDR_q_value <= 0.25,]$sig <- TRUE
sens$Timepoint <- sens$Time
sens$Timepoint <- factor(sens$Timepoint, levels=c("8h","24h"))
sens[sens$FDR_q_value == 0,]$minusLogFDR <- 4
sens[sens$Drug_set == "Bcr.Abl kinase inhibitor",]$Drug_set <- "Bcr-Abl kinase inhibitor"
sens[sens$DrugTreatment=="Trabectidin",]$DrugTreatment <- "Trabectedin"

gsea <- read.csv(synapser::synGet("syn68702141")$path)
gsea$Timepoint <- gsea$Time
gsea$Timepoint <- factor(gsea$Timepoint, levels=c("8h", "24h"))
gsea$individualID <- gsea$MPNST
gsea$Drug <- gsea$DrugTreatment
gsea$signif <- gsea$sig
gsea$feature <- sub("KEGG_","",gsea$Feature_set)
gsea[gsea$FDR_q_value==0,]$minusLogFDR <- 4
keepCols2 <- c("Drug_set", "DrugTreatment","NES","minusLogFDR","sig","Ranking","Timepoint")
#gsea[gsea$DrugTreatment=="Trabectidin",]$DrugTreatment <- "Trabectedin"

tf <- read.csv(synapser::synGet("syn64420968")$path)
tf$feature <- tf$source
tf$minusLogFDR <- -log10(tf$padj)
#tf[tf$Drug=="Trabectidin",]$Drug <- "Trabectedin"

maxAbsNESVals <- list("Differential Expression" = max(abs(diffexp$log2FoldChange), na.rm=TRUE),
                      "Gene Set Enrichment" = max(abs(gsea$NES), na.rm=TRUE),
                      "Transcription Factor Enrichment" = max(abs(tf$score), na.rm=TRUE))
maxLogFDRVals <- list("Differential Expression" = max(diffexp$minusLogFDR, na.rm=TRUE),
                      "Gene Set Enrichment" = max(gsea$minusLogFDR, na.rm=TRUE),
                      "Transcription Factor Enrichment" = max(tf$minusLogFDR, na.rm=TRUE))

#### find drug targets ####
# only interested in MEKi, CDK, HDAC, SHP2
# for each drugTreatment MOA
tested.drugs <- unique(c(diffexp$Drug,gsea$DrugTreatment,sens$DrugTreatment,tf$Drug))
drug.info <- list("Palbociclib" = "CDK inhibitor", "Ribociclib" = "CDK inhibitor",
                  "Trabectedin" = "Chemotherapy", "Ifosfamide" = "DNA alkylating agent",
                  "Decitabine" = "DNMT inhibitor", "Vorinostat" = "HDAC inhibitor",
                  "Mirdametinib" = "MEK inhibitor", "Selumetinib" = "MEK inhibitor",
                  "Trametinib" = "MEK inhibitor", "Capmatinib" = "MET inhibitor",
                  "Olaparib" = "PARP inhibitor", "RMC4630" = "SHP2 inhibitor",
                  "TNO155" = "SHP2 inhibitor", "Doxorubicin" = "TOP inhibitor",
                  #"Irinotecan" = "TOP inhibitor",
                  "Verteporfin" = "YAP inhibitor")

# also have version matching PRISM drug moa sets
drug.info2 <- list("Palbociclib" = "CDK inhibitor", "Ribociclib" = "CDK inhibitor",
                  "Trabectedin" = "Chemotherapy", "Ifosfamide" = "DNA alkylating agent",
                  "Decitabine" = "DNA methyltransferase inhibitor", "Vorinostat" = "HDAC inhibitor",
                  "Mirdametinib" = "MEK inhibitor", "Selumetinib" = "MEK inhibitor",
                  "Trametinib" = "MEK inhibitor", "Capmatinib" = "MET inhibitor",
                  "Olaparib" = "PARP inhibitor", "RMC4630" = "SHP2 inhibitor",
                  "TNO155" = "SHP2 inhibitor", "Doxorubicin" = "Topoisomerase inhibitor",
                  #"Irinotecan" = "Topoisomerase inhibitor",
                  "Verteporfin" = "YAP inhibitor")

all(tested.drugs %in% names(drug.info)) # TRUE

drug.target.info <- read.csv("https://ndownloader.figshare.com/files/20237712")
gen.drug.targets <- list("CDK inhibitor" = "CDK", # default in case drug target info is unavailable
                     "Chemotherapy" = "",
                     "DNA alkylating agent" = "CYP",
                     "DNMT inhibitor" = "DNMT",
                     "HDAC inhibitor" = "HDAC",
                     "MEK inhibitor" = "MAP2K",
                     "MET inhibitor" = "MET",
                     "PARP inhibitor" = "PARP",
                     "SHP2 inhibitor" = "SHP2",
                     "TOP inhibitor" = "TOP",
                     "YAP inhibitor" = c("YAP","TEAD")) # based on https://www.medchemexpress.com/Verteporfin.html?srsltid=AfmBOoqj459z_TwpGATdUezlD2avz3V7tYHYfvbRqodBhC9omcoablXK
drug.targets <- list()
# there is no "TAZ" gene based on grepl
for (i in names(drug.info)) {
  cat(i,"\n")

  # set default targets
  temp.moa <- drug.info[[i]]
  drug.targets[[i]] <- gen.drug.targets[[temp.moa]]

  # try to get specific targets from PRISM
  if (tolower(i) %in% tolower(drug.target.info$name)) {
    cat("in PRISM\n")
    prism.targets <- unique(unlist(strsplit(na.omit(drug.target.info[tolower(drug.target.info$name) == tolower(i),c("name","target")])$target, ", ")[[1]]))
    if (all(!is.na(prism.targets))) {
      cat("targets in PRISM\n")
      drug.targets[[i]] <- prism.targets
    }
  }
}

#### find target effects for each drug ####
measured.targets <- list()
Drug <- names(drug.info)
max.effect <- data.frame(Drug, Gene = NA, Log2FC = NA, p=NA, q=NA, Timepoint = NA, MPNST = NA)
for (i in names(drug.info)){
  temp.targets <- data.frame()
  temp.diffexp <- diffexp[diffexp$hgnc_symbol != "" & diffexp$Drug == i,]
  for (j in drug.targets[[i]]) {
    temp.targets <- rbind(temp.targets, na.omit(temp.diffexp[startsWith(j,temp.diffexp$hgnc_symbol),]))
  }
  if (nrow(temp.targets) > 0) {
    measured.targets[[i]] <- temp.targets
    max.effect[max.effect$Drug == i,]$Timepoint <- levels(temp.targets$Timepoint)[temp.targets[which.max(abs(temp.targets$log2FoldChange)),]$Timepoint]
    max.effect[max.effect$Drug == i,]$Gene <- temp.targets[which.max(abs(temp.targets$log2FoldChange)),]$hgnc_symbol
    max.effect[max.effect$Drug == i,]$Log2FC <- temp.targets[which.max(abs(temp.targets$log2FoldChange)),]$log2FoldChange
    max.effect[max.effect$Drug == i,]$p <- temp.targets[which.max(abs(temp.targets$log2FoldChange)),]$pvalue
    max.effect[max.effect$Drug == i,]$q <- temp.targets[which.max(abs(temp.targets$log2FoldChange)),]$padj
    max.effect[max.effect$Drug == i,]$MPNST <- temp.targets[which.max(abs(temp.targets$log2FoldChange)),]$individualID
  }
}
target.diffexp <- data.table::rbindlist(measured.targets, use.names=FALSE)
write.csv(target.diffexp, "target_diffexp.csv", row.names = FALSE)
write.csv(max.effect, "maximumEffect.csv", row.names=FALSE)
target.diffexp <- read.csv("target_diffexp.csv")
max.effect <- read.csv("maximumEffect.csv")

#### overlapping results across all drugs ####
# DMEA
sig.sens <- na.omit(sens[sens$sig,])
nSig.sens <- plyr::ddply(sig.sens, .(Drug_set), summarize,
                         nSigDrugs = length(unique(DrugTreatment)),
                         nSigMPNST = length(unique(MPNST)))
allSig.sens <- nSig.sens[nSig.sens$nSigDrugs == max(nSig.sens$nSigDrugs) &
                           nSig.sens$nSigMPNST == max(nSig.sens$nSigMPNST),]$Drug_set # 14 drugs in all 3 MPNST
dot.df <- plyr::ddply(sig.sens[sig.sens$Drug_set %in% allSig.sens,], .(Drug_set, Time), summarize,
                      medianNES = median(NES),
                      sdNES = sd(NES))
dot.df$Time <- factor(dot.df$Time, levels=c("8h","24h"))
maxAbsNES <- max(abs(dot.df$medianNES)) # 4.52
maxLogFDR <- max(dot.df$sdNES) # 4
dot.plot <- ggplot2::ggplot(dot.df, ggplot2::aes(x = Time, y = reorder(Drug_set, medianNES),
                                                 color = medianNES, size = sdNES)) +
  ggplot2::geom_point() + ggplot2::scale_x_discrete(limits = c("8h","24h")) +
  scale_color_gradient2(low="blue",high="red", mid="grey",
                        limits=c(-maxAbsNES, maxAbsNES)) +
  theme_classic(base_size=12) + scale_size_continuous(limits=c(0,maxLogFDR)) +
  ggplot2::labs(color = "Median NES", size = "SD NES", title="DMEA") +
  theme(axis.title=element_blank(), plot.title=element_text(hjust=0.5, face="bold"),
        axis.text.x = element_text(angle=45, vjust=1, hjust=1)) +
  geom_point(data = dot.df, col = "black", stroke = 1.5, shape = 21)
customHeight <- ifelse(length(unique(dot.df$Drug_set))/10 < 2, 2, length(unique(dot.df$Drug_set))/10)
ggplot2::ggsave("DMEA_dotPlot.pdf", dot.plot, width=4, height=customHeight)

dot.plot <- ggplot2::ggplot(dot.df, ggplot2::aes(x = Time, y = reorder(Drug_set, -medianNES),
                                                 color = medianNES, size = sdNES)) +
  ggplot2::geom_point() + ggplot2::scale_x_discrete(limits = c("8h","24h")) +
  scale_color_gradient2(low="blue",high="red", mid="grey",
                        limits=c(-maxAbsNES, maxAbsNES)) +
  theme_classic(base_size=12) + scale_size_continuous(limits=c(0,maxLogFDR)) +
  ggplot2::labs(color = "Median NES", size = "SD NES", title="DMEA") +
  theme(axis.title=element_blank(), plot.title=element_text(hjust=0.5, face="bold"),
        axis.text.x = element_text(angle=45, vjust=1, hjust=1)) +
  geom_point(data = dot.df, col = "black", stroke = 1.5, shape = 21)
customHeight <- ifelse(length(unique(dot.df$Drug_set))/10 < 2, 2, length(unique(dot.df$Drug_set))/10)
ggplot2::ggsave("DMEA_dotPlot_reverseOrder.pdf", dot.plot, width=4, height=customHeight)

# diffexp
dot.df <- na.omit(diffexp[diffexp$signif,])
rep.genes <- na.omit(unique(dot.df$hgnc_symbol[duplicated(dot.df$hgnc_symbol)])) # 15527 genes - too many to plot so just focus on targets
mean.dot.df <- plyr::ddply(dot.df, .(hgnc_symbol), summarize,
                           N_drugs = length(na.omit(unique(Drug))),
                           drugs = paste0(na.omit(unique(Drug),), collapse=", "),
                           N_MPNST = length(na.omit(unique(individualID))),
                           MPNST = paste0(na.omit(unique(individualID),), collapse=", "))
mean.dot.df$all_same <- mean.dot.df$all_dn | mean.dot.df$all_up
rep.genes <- na.omit(unique(mean.dot.df[mean.dot.df$N_drugs==max(mean.dot.df$N_drugs) &
                                          mean.dot.df$N_MPNST==max(mean.dot.df$N_MPNST),]$hgnc_symbol)) # DUSP6; not all the same dir
dot.df <- diffexp[diffexp$hgnc_symbol %in% rep.genes,]
dot.df$Time <- dot.df$Timepoint
dot.df$Drug_set <- dot.df$hgnc_symbol
dot.df$MPNST <- dot.df$individualID
dot.df$Time <- factor(dot.df$Time, levels=c("8h","24h"))
dot.df$sig <- dot.df$signif
dot.df$DrugTreatment <- dot.df$Drug
dot.df$NES <- dot.df$log2FoldChange
maxAbsNES <- ceiling(max(abs(na.omit(dot.df$NES)))) # 5
maxAbsP <- ceiling(max(na.omit(dot.df$minusLogFDR))) # 209
DMEAdotPlots(dot.df,
             fname = paste0("DUSP6","_diffexp_dotPlot"),
             maxAbsNES = maxAbsNES, maxLogFDR = maxAbsP,
             colorLab="Log2FC", width=12, height=3)

# TF
dot.df <- tf[tf$signif,]
rep.genes <- na.omit(unique(dot.df$source[duplicated(dot.df$source)])) # 103
mean.dot.df <- plyr::ddply(dot.df, .(source), summarize,
                           N_drugs = length(na.omit(unique(Drug))),
                           drugs = paste0(na.omit(unique(Drug),), collapse=", "),
                           N_MPNST = length(na.omit(unique(individualID))),
                           MPNST = paste0(na.omit(unique(individualID),), collapse=", "),
                           all_up = all(score > 0),
                           all_dn = all(score < 0))
mean.dot.df$all_same <- mean.dot.df$all_dn | mean.dot.df$all_up
rep.genes <- na.omit(unique(mean.dot.df[mean.dot.df$N_drugs==max(mean.dot.df$N_drugs) &
                                          mean.dot.df$N_MPNST==max(mean.dot.df$N_MPNST),]$source)) # SP1; not all same dir
dot.df <- tf[tf$source %in% rep.genes,]
dot.df$Time <- dot.df$Timepoint
dot.df$Drug_set <- dot.df$source
dot.df$MPNST <- dot.df$individualID
dot.df$Time <- factor(dot.df$Time, levels=c("8h","24h"))
dot.df$sig <- dot.df$signif
dot.df$DrugTreatment <- dot.df$Drug
dot.df$NES <- dot.df$score
dot.df$minusLogFDR <- -log10(dot.df$padj)
maxAbsNES <- ceiling(max(abs(na.omit(dot.df$NES)))) # 16
maxAbsP <- ceiling(max(na.omit(dot.df$minusLogFDR))) # 51
DMEAdotPlots(dot.df,
             fname = paste0("SP1","_TF_dotPlot"),
             maxAbsNES = maxAbsNES, maxLogFDR = maxAbsP,
             width=12, height=3,
             colorLab="Score")

# GSEA
dot.df <- gsea[gsea$sig,]
rep.genes <- na.omit(unique(dot.df$feature[duplicated(dot.df$feature)])) # 49
mean.dot.df <- plyr::ddply(dot.df, .(feature), summarize,
                           N_drugs = length(na.omit(unique(DrugTreatment))),
                           drugs = paste0(na.omit(unique(DrugTreatment),), collapse=", "),
                           N_MPNST = length(na.omit(unique(MPNST))),
                           MPNST = paste0(na.omit(unique(MPNST),), collapse=", "),
                           all_up = all(NES > 0),
                           all_dn = all(NES < 0))
mean.dot.df$all_same <- mean.dot.df$all_dn | mean.dot.df$all_up
rep.genes <- na.omit(unique(mean.dot.df[mean.dot.df$N_drugs==max(mean.dot.df$N_drugs) &
                                          mean.dot.df$N_MPNST==max(mean.dot.df$N_MPNST),]$feature)) # 5

dot.df <- gsea[gsea$feature %in% rep.genes,]
dot.df$Drug_set <- sub("HALLMARK_","",dot.df$feature)
dot.df$Time <- factor(dot.df$Time, levels=c("8h","24h"))
maxAbsNES <- max(abs(dot.df$NES)) # 3.26
maxLogFDR <- max(dot.df$minusLogFDR) # 4
DMEAdotPlots(dot.df,
             fname = paste0("GSEA_dotPlot"),
             maxAbsNES = maxAbsNES, maxLogFDR = maxLogFDR,
             width=12, height=4)
