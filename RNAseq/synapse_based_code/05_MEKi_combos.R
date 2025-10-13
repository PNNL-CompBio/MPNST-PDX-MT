# compile DMEA for potential synergy
rm(list=ls())
library(plyr);library(dplyr);library(synapser);library(data.table);library(ggplot2);library(patchwork)
synapser::synLogin()
setwd("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/MPNST-PDX-MT/RNAseq/synapse_based_code")
dir.create("comboPrediction_MEK")
setwd("comboPrediction_MEK")

DMEAdotPlots <- function(dot.df,  title="", fname, min.width = 5, min.height = 4, 
                         maxAbsNES = max(abs(dot.df[,color])), 
                         maxLogFDR = max(dot.df[,size]), x.order=sort(unique(dot.df$DrugTreatment)),
                         colorLab = "NES") {
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
  customHeight <- ifelse(length(unique(dot.df$Drug_set))/5 < min.height, min.height, length(unique(dot.df$Drug_set))/5)
  ggplot2::ggsave(paste0(fname,".pdf"), dot.plot, width=min.width, height=customHeight)
  ggplot2::ggsave(paste0(fname,"_taller.pdf"), dot.plot, width=min.width, height=(customHeight+1))
  ggplot2::ggsave(paste0(fname,"_widerTaller.pdf"), dot.plot, width=(min.width+2), height=(customHeight+2))
  ggplot2::ggsave(paste0(fname,"_evenWiderTaller.pdf"), dot.plot, width=(min.width+3), height=(customHeight+3))
  ggplot2::ggsave(paste0(fname,"_evenTaller.pdf"), dot.plot, width=min.width, height=(customHeight+3))
  ggplot2::ggsave(paste0(fname,"_evenWider.pdf"), dot.plot, width=(min.width+3), height=(customHeight+2))
  
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
diffexp <- read.csv(synapser::synGet("syn69910880")$path) # read.csv(synapser::synGet('syn64397743')$path)
diffexp <- read.csv(synapser::synGet('syn64397743',version=19)$path)
write.csv(diffexp,"totalRNASeqDiffExtoDate.csv",row.names=FALSE)
diffexp$log2FoldChange <- as.numeric(diffexp$log2FoldChange)
diffexp$padj <- as.numeric(diffexp$padj)
diffexp$Timepoint <- diffexp$Time
diffexp$Timepoint <- factor(diffexp$Timepoint, levels=c("8h", "24h"))
diffexp$feature <- diffexp$hgnc_symbol
diffexp$minusLogFDR <- -log10(diffexp$padj)
diffexp <- diffexp[!is.na(diffexp$hgnc_symbol) & diffexp$hgnc_symbol != "",]
diffexp[diffexp$Drug=="Trabectidin",]$Drug <- "Trabectedin"

sens <- read.csv(synapser::synGet("syn69910844")$path) # read.csv(synapser::synGet("syn65672258")$path)
sens$sig <- FALSE
sens[sens$p_value <= 0.05 & sens$FDR_q_value <= 0.25,]$sig <- TRUE
sens$Timepoint <- sens$Time
sens$Timepoint <- factor(sens$Timepoint, levels=c("8h","24h"))
sens[sens$FDR_q_value == 0,]$minusLogFDR <- 4
sens[sens$Drug_set == "Bcr.Abl kinase inhibitor",]$Drug_set <- "Bcr-Abl kinase inhibitor"
sens[sens$DrugTreatment=="Trabectidin",]$DrugTreatment <- "Trabectedin"

gsea <- read.csv(synapser::synGet("syn69910859")$path) # read.csv(synapser::synGet("syn68702141")$path)
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
doi <- drug.info[grepl("MEK", drug.info) | grepl("CDK", drug.info) |
                   grepl("TOP", drug.info)| grepl("HDAC", drug.info) | grepl("SHP2", drug.info)]
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
doi2 <- drug.info2[grepl("MEK", drug.info2) | grepl("CDK", drug.info2) |
                    grepl("Topoisomerase", drug.info2) |
                    grepl("HDAC", drug.info2)]

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

#### overlapping results across MEK, CDK, HDAC, TOP, SHP2 ####
# DMEA
dot.df <- na.omit(sens[sens$DrugTreatment %in% names(doi2) & sens$Drug_set %in% doi2,])
dot.df$Time <- factor(dot.df$Time, levels=c("8h","24h"))
maxAbsNES <- max(abs(dot.df$NES)) # 4.52
maxLogFDR <- max(dot.df$minusLogFDR) # 4
DMEAdotPlots(dot.df, 
             fname = paste0("MEK-HDAC-TOP-CDK","_DMEA_dotPlot"), 
             maxAbsNES = maxAbsNES, maxLogFDR = maxLogFDR, 
             x.order=c(names(doi[doi=="MEK inhibitor"]), 
                       names(doi[doi=="HDAC inhibitor"]),
                       names(doi[doi=="TOP inhibitor"]),
                       names(doi[doi=="CDK inhibitor"])))

# diffexp
dot.df <- na.omit(diffexp[diffexp$Drug %in% names(doi2) & diffexp$signif,])
dot.df$moa <- ""
moi <- paste(c("MEK","HDAC","TOP","CDK"),"inhibitor")
for (i in moi) {
  dot.df[dot.df$Drug %in% names(doi[doi==i]),]$moa <- i
}
rep.genes <- na.omit(unique(dot.df$hgnc_symbol[duplicated(dot.df$hgnc_symbol)])) # 4149 genes - too many to plot so just focus on targets
mean.dot.df <- plyr::ddply(dot.df, .(hgnc_symbol, Timepoint, individualID), summarize,
                           N_moa_up = length(na.omit(unique(moa[log2FoldChange>0]))),
                           moas_up = paste0(na.omit(unique(moa[log2FoldChange>0])), collapse=", "),
                           N_moa_dn = length(na.omit(unique(moa[log2FoldChange<0]))),
                           moas_dn = paste0(na.omit(unique(moa[log2FoldChange<0])), collapse=", "),
                           N_moa = length(na.omit(unique(moa))),
                           moas = paste0(na.omit(unique(moa),), collapse=", "),
                           all_up = all(log2FoldChange > 0),
                           all_dn = all(log2FoldChange < 0))
mean.dot.df$all_same <- mean.dot.df$all_dn | mean.dot.df$all_up
rep.genes <- na.omit(unique(mean.dot.df[mean.dot.df$N_moa>1,]$hgnc_symbol)) # 1608 - too many
rep.genes <- na.omit(unique(mean.dot.df[mean.dot.df$N_moa>1 & mean.dot.df$all_same,]$hgnc_symbol)) # 890 - too many
rep.genes <- na.omit(unique(mean.dot.df[mean.dot.df$N_moa>2,]$hgnc_symbol)) # 131 - too many
rep.genes <- na.omit(unique(mean.dot.df[mean.dot.df$N_moa>2 & mean.dot.df$all_same,]$hgnc_symbol)) # 60 - too many
rep.genes <- na.omit(unique(mean.dot.df[mean.dot.df$N_moa>3,]$hgnc_symbol)) # 6
rep.genes <- na.omit(unique(mean.dot.df[mean.dot.df$N_moa>3 & mean.dot.df$all_same,]$hgnc_symbol)) # 2 too few
rep.genes <- na.omit(unique(mean.dot.df[mean.dot.df$N_moa>3,]$hgnc_symbol)) # 6
dot.df <- diffexp[diffexp$Drug %in% names(doi2) & diffexp$hgnc_symbol %in% rep.genes,]
dot.df$Time <- dot.df$Timepoint
dot.df$Drug_set <- dot.df$hgnc_symbol
dot.df$MPNST <- dot.df$individualID
dot.df$Time <- factor(dot.df$Time, levels=c("8h","24h"))
dot.df$sig <- dot.df$signif
dot.df$DrugTreatment <- dot.df$Drug
dot.df$NES <- dot.df$log2FoldChange
maxAbsNES <- ceiling(max(abs(na.omit(dot.df$NES)))) # 8
maxAbsP <- ceiling(max(na.omit(dot.df$minusLogFDR))) # 209
DMEAdotPlots(dot.df, 
             fname = paste0("MEK-HDAC-TOP-CDK","_diffexpMin4SigMOAs_dotPlot"), 
             maxAbsNES = maxAbsNES, maxLogFDR = maxAbsP, 
             x.order=c(names(doi[doi=="MEK inhibitor"]), 
                       names(doi[doi=="HDAC inhibitor"]),
                       names(doi[doi=="TOP inhibitor"]),
                       names(doi[doi=="CDK inhibitor"]),
                       names(doi[doi=="SHP2 inhibitor"])),
             colorLab="Log2FC")

# target diffexp
target.genes.of.interest <- target.diffexp[target.diffexp$Drug %in% names(doi2),]$hgnc_symbol
dot.df <- diffexp[diffexp$Drug %in% names(doi2) & diffexp$hgnc_symbol %in% target.genes.of.interest,]
dot.df$Time <- dot.df$Timepoint
dot.df$Drug_set <- dot.df$hgnc_symbol
dot.df$MPNST <- dot.df$individualID
dot.df$Time <- factor(dot.df$Time, levels=c("8h","24h"))
dot.df$sig <- dot.df$signif
dot.df$DrugTreatment <- dot.df$Drug
dot.df$NES <- dot.df$log2FoldChange
maxAbsNES <- ceiling(max(abs(na.omit(dot.df$NES)))) # 16
maxAbsP <- ceiling(max(na.omit(dot.df$minusLogFDR))) # 22
DMEAdotPlots(dot.df, 
             fname = paste0("MEK-HDAC-TOP-CDK","_targetDiffexp_dotPlot"), 
             maxAbsNES = maxAbsNES, maxLogFDR = maxAbsP, 
             x.order=c(names(doi[doi=="MEK inhibitor"]), 
                       names(doi[doi=="HDAC inhibitor"]),
                       names(doi[doi=="TOP inhibitor"]),
                       names(doi[doi=="CDK inhibitor"])),
             colorLab="Log2FC")

# target diffexp just based on gene symbols
dot.df <- diffexp[diffexp$Drug %in% names(doi2) & 
                    (startsWith(diffexp$hgnc_symbol,"MAP2K") |
                       startsWith(diffexp$hgnc_symbol, "HDAC") |
                       startsWith(diffexp$hgnc_symbol,"TOP") |
                       startsWith(diffexp$hgnc_symbol,"CDK")) & diffexp$signif,]
dot.df <- diffexp[diffexp$Drug %in% names(doi2) & diffexp$hgnc_symbol %in% dot.df$hgnc_symbol,]
dot.df$Time <- dot.df$Timepoint
dot.df$Drug_set <- dot.df$hgnc_symbol
dot.df$MPNST <- dot.df$individualID
dot.df$Time <- factor(dot.df$Time, levels=c("8h","24h"))
dot.df$sig <- dot.df$signif
dot.df$DrugTreatment <- dot.df$Drug
dot.df$NES <- dot.df$log2FoldChange
maxAbsNES <- ceiling(max(abs(na.omit(dot.df$NES)))) # 16
maxAbsP <- ceiling(max(na.omit(dot.df$minusLogFDR))) # 22
DMEAdotPlots(dot.df, 
             fname = paste0("MEK-HDAC-TOP-CDK","_targetDiffexp_dotPlot_v2"), 
             maxAbsNES = maxAbsNES, maxLogFDR = maxAbsP, 
             x.order=c(names(doi[doi=="MEK inhibitor"]), 
                       names(doi[doi=="HDAC inhibitor"]),
                       names(doi[doi=="TOP inhibitor"]),
                       names(doi[doi=="CDK inhibitor"])),
             colorLab="Log2FC")

expected.targets <- c(paste0("MAP2K",seq(1,19)),paste0("CDK",seq(1,19)),
                      paste0("HDAC",seq(1,19)),paste0("TOP",seq(1,19)))
dot.df <- na.omit(diffexp[diffexp$Drug %in% names(doi2) & 
                    diffexp$hgnc_symbol %in% expected.targets & diffexp$signif,])
dot.df <- diffexp[diffexp$Drug %in% names(doi2) & 
                    diffexp$hgnc_symbol %in% dot.df$hgnc_symbol,]
dot.df$Time <- dot.df$Timepoint
dot.df$Drug_set <- dot.df$hgnc_symbol
dot.df$MPNST <- dot.df$individualID
dot.df$Time <- factor(dot.df$Time, levels=c("8h","24h"))
dot.df$sig <- dot.df$signif
dot.df$DrugTreatment <- dot.df$Drug
dot.df$NES <- dot.df$log2FoldChange
maxAbsNES <- ceiling(max(abs(na.omit(dot.df$NES)))) # 16
maxAbsP <- ceiling(max(na.omit(dot.df$minusLogFDR))) # 22
DMEAdotPlots(dot.df, 
             fname = paste0("MEK-HDAC-TOP-CDK","_targetDiffexp_dotPlot_v3"), 
             maxAbsNES = maxAbsNES, maxLogFDR = maxAbsP, 
             x.order=c(names(doi[doi=="MEK inhibitor"]), 
                       names(doi[doi=="HDAC inhibitor"]),
                       names(doi[doi=="TOP inhibitor"]),
                       names(doi[doi=="CDK inhibitor"])),
             colorLab="Log2FC")

# just mirda and vorinostat
expected.targets <- c(paste0("MAP2K",seq(1,19)),
                      paste0("HDAC",seq(1,19)))
dot.df <- na.omit(diffexp[diffexp$Drug %in% c("Mirdametinib","Vorinostat") & 
                            diffexp$hgnc_symbol %in% expected.targets & diffexp$signif,])
dot.df <- diffexp[diffexp$Drug %in% c("Mirdametinib","Vorinostat") & 
                    diffexp$hgnc_symbol %in% dot.df$hgnc_symbol,]
dot.df$Time <- dot.df$Timepoint
dot.df$Drug_set <- dot.df$hgnc_symbol
dot.df$MPNST <- dot.df$individualID
dot.df$Time <- factor(dot.df$Time, levels=c("8h","24h"))
dot.df$sig <- dot.df$signif
dot.df$DrugTreatment <- dot.df$Drug
dot.df$NES <- dot.df$log2FoldChange
maxAbsNES <- ceiling(max(abs(na.omit(dot.df$NES)))) # 16
maxAbsP <- ceiling(max(na.omit(dot.df$minusLogFDR))) # 22
DMEAdotPlots(dot.df, 
             fname = paste0("mirdaVorin","_targetDiffexp_dotPlot_v3"), 
             maxAbsNES = maxAbsNES, maxLogFDR = maxAbsP, 
             x.order=c("Mirdametinib","Vorinostat"),
             colorLab="Log2FC")

# TF
dot.df <- tf[tf$Drug %in% names(doi2) & tf$signif,]
dot.df$moa <- ""
moi <- paste(c("MEK","HDAC","TOP","CDK"),"inhibitor")
for (i in moi) {
  dot.df[dot.df$Drug %in% names(doi[doi==i]),]$moa <- i
}
rep.genes <- na.omit(unique(dot.df$source[duplicated(dot.df$source)])) # 86
mean.dot.df <- plyr::ddply(dot.df, .(source, Timepoint, individualID), summarize,
                           N_moa_up = length(na.omit(unique(moa[score>0]))),
                           moas_up = paste0(na.omit(unique(moa[score>0])), collapse=", "),
                           N_moa_dn = length(na.omit(unique(moa[score<0]))),
                           moas_dn = paste0(na.omit(unique(moa[score<0])), collapse=", "),
                           N_moa = length(na.omit(unique(moa))),
                           moas = paste0(na.omit(unique(moa),), collapse=", "),
                           all_up = all(score > 0),
                           all_dn = all(score < 0))
mean.dot.df$all_same <- mean.dot.df$all_dn | mean.dot.df$all_up
rep.genes <- na.omit(unique(mean.dot.df[mean.dot.df$N_moa>1,]$source)) # 56
rep.genes <- na.omit(unique(mean.dot.df[mean.dot.df$N_moa>1 & mean.dot.df$all_same,]$source)) # 49
rep.genes <- na.omit(unique(mean.dot.df[mean.dot.df$N_moa>2,]$source)) # 25
rep.genes <- na.omit(unique(mean.dot.df[mean.dot.df$N_moa>2 & mean.dot.df$all_same,]$source)) # 25
rep.genes <- na.omit(unique(mean.dot.df[mean.dot.df$N_moa>3,]$source)) # 0
rep.genes <- na.omit(unique(mean.dot.df[mean.dot.df$N_moa>2,]$source)) # 25
dot.df <- tf[tf$Drug %in% names(doi2) & tf$source %in% rep.genes,]
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
             fname = paste0("MEK-HDAC-CDK","_TFmin3SigMOAs_dotPlot"), 
             maxAbsNES = maxAbsNES, maxLogFDR = maxAbsP, 
             x.order=c(names(doi[doi=="MEK inhibitor"]), 
                       names(doi[doi=="HDAC inhibitor"]),
                       names(doi[doi=="TOP inhibitor"]),
                       names(doi[doi=="CDK inhibitor"]),
                       names(doi[doi=="SHP2 inhibitor"])),
             colorLab="Score")

# GSEA
dot.df <- gsea[gsea$DrugTreatment %in% names(doi2) & gsea$sig,]
dot.df$moa <- ""
moi <- paste(c("MEK","HDAC","TOP","CDK"),"inhibitor")
for (i in moi) {
  dot.df[dot.df$DrugTreatment %in% names(doi[doi==i]),]$moa <- i
}
rep.genes <- na.omit(unique(dot.df$feature[duplicated(dot.df$feature)])) # 38
mean.dot.df <- plyr::ddply(dot.df, .(feature, Time, MPNST), summarize,
                           N_moa_up = length(na.omit(unique(moa[NES>0]))),
                           moas_up = paste0(na.omit(unique(moa[NES>0])), collapse=", "),
                           N_moa_dn = length(na.omit(unique(moa[NES<0]))),
                           moas_dn = paste0(na.omit(unique(moa[NES<0])), collapse=", "),
                           N_moa = length(na.omit(unique(moa))),
                           moas = paste0(na.omit(unique(moa),), collapse=", "),
                           all_up = all(NES > 0),
                           all_dn = all(NES < 0))
mean.dot.df$all_same <- mean.dot.df$all_dn | mean.dot.df$all_up
rep.genes <- na.omit(unique(mean.dot.df[mean.dot.df$N_moa>1,]$feature)) # 28
rep.genes <- na.omit(unique(mean.dot.df[mean.dot.df$N_moa>1 & mean.dot.df$all_same,]$feature)) # 23
rep.genes <- na.omit(unique(mean.dot.df[mean.dot.df$N_moa>2,]$feature)) # 11
rep.genes <- na.omit(unique(mean.dot.df[mean.dot.df$N_moa>2 & mean.dot.df$all_same,]$feature)) # 9
rep.genes <- na.omit(unique(mean.dot.df[mean.dot.df$N_moa>3,]$feature)) # 2
rep.genes <- na.omit(unique(mean.dot.df[mean.dot.df$N_moa>3 & mean.dot.df$all_same,]$feature)) # 0
rep.genes <- na.omit(unique(mean.dot.df[mean.dot.df$N_moa>2,]$feature)) # 11

dot.df <- gsea[gsea$DrugTreatment %in% names(doi2) & gsea$feature %in% rep.genes,]
dot.df$Drug_set <- dot.df$feature
dot.df$Time <- factor(dot.df$Time, levels=c("8h","24h"))
maxAbsNES <- max(abs(dot.df$NES)) # 4.52; 8
maxLogFDR <- max(dot.df$minusLogFDR) # 4; 4
DMEAdotPlots(dot.df, 
             fname = paste0("MEK-HDAC-TOP-CDK","_GSEAmin3SigMOAs_dotPlot"), 
             maxAbsNES = maxAbsNES, maxLogFDR = maxLogFDR, 
             x.order=c(names(doi[doi=="MEK inhibitor"]), 
                       names(doi[doi=="HDAC inhibitor"]),
                       names(doi[doi=="TOP inhibitor"]),
                       names(doi[doi=="CDK inhibitor"]),
                       names(doi[doi=="SHP2 inhibitor"])))
# 
# #### PDX data: MEKi+HDACi ####
# #tumor.size <- read.csv("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/MPNST-PDX-MT/Wu225 Mirda + Vorinostat.csv")
# tumor.size <- read.csv(synapser::synGet("syn68900596")$path)
# 
# # get colors for mirda and vorinostat
# library(RColorBrewer); library(scales)
# cols1=c("#000000",brewer.pal(8,'Dark2'),brewer.pal(15-8,'Set2'),"mediumorchid1", "cornflowerblue", "#004B4B", "#4B0026")
# drug.info <- list("palbociclib" = "CDK inhibitor", "ribociclib" = "CDK inhibitor",
#                   "trabectedin" = "Chemotherapy", "Ifosfamide" = "DNA alkylating agent",
#                   "decitabine" = "DNMT inhibitor", "vorinostat" = "HDAC inhibitor",
#                   "mirdametinib" = "MEK inhibitor", "selumetinib" = "MEK inhibitor",
#                   "trametinib" = "MEK inhibitor", "capmatinib" = "MET inhibitor",
#                   "olaparib" = "PARP inhibitor", "RMC4630" = "SHP2 inhibitor",
#                   "TNO155" = "SHP2 inhibitor", "doxorubicin" = "TOP inhibitor",
#                   "irinotecan" = "TOP inhibitor", "verteporfin" = "YAP inhibitor",
#                   "pexidartinib" = "KIT inhibitor", "IAG933" = "YAP inhibitor", "SN38" = "TOP inhibitor")
# names(cols1) <- c("DMSO", names(drug.info))
# scales::show_col(cols1[c("DMSO","mirdametinib","vorinostat")]) # black, brown, yellow?
# ggplot(tumor.size,aes(x=`Treatment.Day..`,y=`Mean.Size`, color=Treatment)) + 
#   geom_point() + geom_errorbar(aes(ymin=`Mean.Size`-SEM, ymax=`Mean.Size`+SEM))+
#   geom_smooth(se=FALSE, linetype="dashed")+
#   theme_classic() + labs(y="Mean Tumor Size", x="Treatment Duration (Days)") +
#   scale_color_manual(values=c("black","red","blue","forestgreen"), breaks=c("Vehicle","Mirda","Vorinostat","Mirda + Vorinostat"))
# setwd("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/MPNST-PDX-MT")
# ggsave("WU225_mirdaVorinostat_PDX.pdf", width=6, height=3) # was width 4
# 
# combo.t <- t.test(tumor.size[tumor.size$Treatment %in% c("Mirda","Vorinostat"),]$`Mean.Size`,
#                   tumor.size[tumor.size$Treatment == "Mirda + Vorinostat",]$`Mean.Size`, "greater")
# combo.t # p = 0.000649
# 
# mek.combo.t <- t.test(tumor.size[tumor.size$Treatment == "Mirda",]$`Mean.Size`,
#                   tumor.size[tumor.size$Treatment == "Mirda + Vorinostat",]$`Mean.Size`, "greater", paired=TRUE)
# mek.combo.t # p = 0.001507
# 
# vor.combo.t <- t.test(tumor.size[tumor.size$Treatment == "Vorinostat",]$`Mean.Size`,
#                       tumor.size[tumor.size$Treatment == "Mirda + Vorinostat",]$`Mean.Size`, "greater", paired=TRUE)
# vor.combo.t # p = 0.002504
# 
# veh.mek.t <- t.test(tumor.size[tumor.size$Treatment == "Vehicle",]$`Mean.Size`,
#                       tumor.size[tumor.size$Treatment == "Mirda",]$`Mean.Size`, "greater", paired=TRUE)
# veh.mek.t # p = 0.006131
# 
# veh.vor.t <- t.test(tumor.size[tumor.size$Treatment == "Vehicle",]$`Mean.Size`,
#                       tumor.size[tumor.size$Treatment == "Vorinostat",]$`Mean.Size`, "greater", paired=TRUE)
# veh.vor.t # p = 0.004812
# 
# #### mirda combos predicted ####
# mirda.sens <- sens[sens$DrugTreatment=="Mirdametinib" & sens$sig & sens$NES < 0,]
# mirda.mean <- plyr::ddply(mirda.sens, .(Drug_set), summarize,
#                           meanNES=mean(NES),
#                           medianNES=median(NES),
#                           sdNES=sd(NES),
#                           N=dplyr::n(),
#                           N_MPNST = length(unique(MPNST)),
#                           MPNST = paste0(sort(unique(MPNST)), collapse=", "))
# mirda.mean$Tested <- "plain"
# mirda.mean[mirda.mean$Drug_set %in% drug.info2,]$Tested <- "bold"
# ggplot(mirda.mean, aes(y=reorder(Drug_set, -meanNES), x=meanNES, fill=MPNST)) + 
#   geom_bar(stat="identity") + theme_classic() + scale_y_discrete(position="right") +
#   theme(axis.title.y=element_blank(), legend.position="bottom") + labs(x="Mean NES") +
#   scale_fill_manual(values=RColorBrewer::brewer.pal(length(unique(mirda.mean$MPNST)),"Dark2"), breaks=unique(mirda.mean$MPNST))
#   #scale_fill_manual(values=c("black","darkgrey","lightgrey"), breaks=c("JH-2-002 & MN-2","JH-2-002","MN-2"))
# ggsave("Mirda_DMEA_predictions.pdf", width=5, height=5)
# 
# ggplot(mirda.mean, aes(x=reorder(Drug_set, meanNES), y=meanNES, fill=MPNST)) + 
#   geom_bar(stat="identity") + theme_classic() + 
#   theme(axis.title.x=element_blank(), legend.position="top", axis.text.x=element_text(angle=45, vjust=1, hjust=1)) + labs(y="Mean NES") +
#   scale_fill_manual(values=RColorBrewer::brewer.pal(length(unique(mirda.mean$MPNST)),"Dark2"), breaks=unique(mirda.mean$MPNST))
# ggsave("Mirda_DMEA_predictions_horizontal.pdf", width=5, height=5)
# 
# mirda.mean[mirda.mean$Drug_set == "protein tyrosine kinase inhibitor",]$Drug_set <- "Protein TKI"
# mirda.mean[mirda.mean$Drug_set == "tyrosine kinase inhibitor",]$Drug_set <- "TKI"
# mirda.mean[grepl("PDGFR", mirda.mean$Drug_set),]$Drug_set <- "PDGFRi"
# mirda.mean[grepl("Aurora", mirda.mean$Drug_set),]$Drug_set <- "AURKi"
# mirda.mean[grepl("bromodomain", mirda.mean$Drug_set),]$Drug_set <- "BRDi"
# mirda.mean[grepl("adenosine", mirda.mean$Drug_set),]$Drug_set <- "P1i"
# mirda.mean[grepl("glutamate", mirda.mean$Drug_set),]$Drug_set <- "GluRi"
# mirda.mean[grepl("src", mirda.mean$Drug_set),]$Drug_set <- "SRCi"
# mirda.mean[grepl("retinoid ", mirda.mean$Drug_set),]$Drug_set <- "RR Agonist"
# mirda.mean[grepl("inhibitor", mirda.mean$Drug_set),]$Drug_set <- sub(" inhibitor","i",mirda.mean[grepl("inhibitor", mirda.mean$Drug_set),]$Drug_set)
# ggplot(mirda.mean, aes(x=reorder(Drug_set, meanNES), y=meanNES, fill=MPNST)) + 
#   geom_bar(stat="identity") + theme_classic() + 
#   theme(axis.title.x=element_blank(), legend.position="top", axis.text.x=element_text(angle=45, vjust=1, hjust=1)) + labs(y="Mean NES") +
#   scale_fill_manual(values=RColorBrewer::brewer.pal(length(unique(mirda.mean$MPNST)),"Dark2"), breaks=unique(mirda.mean$MPNST))
# ggsave("Mirda_DMEA_predictions_horizontal_shortNames.pdf", width=5, height=3)
# ggsave("Mirda_DMEA_predictions_horizontal_shortNames_lessWide.pdf", width=4, height=3)
# 
# ggplot(mirda.mean, aes(y=reorder(Drug_set, -meanNES), x=meanNES, fill=MPNST)) + 
#   geom_bar(stat="identity") + theme_classic() + scale_y_discrete(position="right") +
#   theme(axis.title.y=element_blank(), legend.position="bottom") + labs(x="Mean NES") +
#   scale_fill_manual(values=RColorBrewer::brewer.pal(length(unique(mirda.mean$MPNST)),"Dark2"), breaks=unique(mirda.mean$MPNST))
# ggsave("Mirda_DMEA_predictions_shortNames.pdf", width=4, height=5)
# ggsave("Mirda_DMEA_predictions_shortNames_smaller.pdf", width=3, height=4)
