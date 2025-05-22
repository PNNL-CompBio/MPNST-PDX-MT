# compile DMEA for potential synergy
rm(list=ls())
library(plyr);library(dplyr);library(synapser);library(data.table);library(ggplot2);library(patchwork)
synapser::synLogin()
setwd("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/MPNST-PDX-MT/RNAseq/synapse_based_code")
dir.create("comboPrediction")
setwd("comboPrediction")

DMEAdotPlots <- function(dot.df,  y.order, title, fname, min.width = 5, min.height = 4, 
                         maxAbsNES = max(abs(dot.df[,color])), 
                         maxLogFDR = max(dot.df[,size])) {
  dot.plot <- ggplot2::ggplot(dot.df, ggplot2::aes(x = DrugTreatment, y = Drug_set, 
                                                   color = NES, size = minusLogFDR)) + 
    facet_grid(Timepoint ~ Ranking) + ggplot2::geom_point() +
    ggplot2::scale_y_discrete(limits = y.order) +
    scale_color_gradient2(low="blue",high="red", mid="grey", 
                          limits=c(-maxAbsNES, maxAbsNES)) +
    theme_classic(base_size=12) + scale_size_continuous(limits=c(0,maxLogFDR)) + 
    ggplot2::labs(color = "NES", size = "-log(FDR)", title=title) +  
    theme(axis.title=element_blank(), plot.title=element_text(hjust=0.5, face="bold"),
          axis.text.x = element_text(angle=45, vjust=1, hjust=1)) +
    geom_point(data = subset(dot.df, sig), col = "black", stroke = 1.5, shape = 21)
  customHeight <- ifelse(length(y.order)/5 < min.height, min.height, length(y.order)/5)
  ggplot2::ggsave(paste0(fname,".pdf"), dot.plot, width=min.width, height=customHeight)
  
  if (all(grepl(" inhibitor", y.order))) {
    dot.df$Inhibitor <- sub(" inhibitor", "", dot.df$Drug_set)
    sigOrder2 <- sub(" inhibitor", "", y.order)
    longNames <- c("Aurora kinase", "Bromodomain", "Sterol demethylase", 
                   "Topoisomerase", "Tyrosine kinase", "Tubulin",
                   "Ribonucleotide reductase")
    names(longNames) <- c("AURK", "BET", "CYP51", "TOP", "TK", "TUB", "RNR")
    if (any(dot.df$Inhibitor %in% longNames)) {
      for (tempName in names(longNames)) {
        if (any(dot.df$Inhibitor == longNames[[tempName]])) {
          dot.df[dot.df$Inhibitor == longNames[[tempName]],]$Inhibitor <- tempName
          sigOrder2 <- sub(longNames[[tempName]], tempName, sigOrder2)
        }
      }
    }
    dot.plot <- ggplot2::ggplot(
      dot.df,
      ggplot2::aes(
        x = DrugTreatment, y = Inhibitor, color = NES,
        size = minusLogFDR
      )
    ) + facet_grid(Timepoint ~ Ranking) +
      ggplot2::geom_point() +
      ggplot2::scale_y_discrete(limits = sigOrder2) +
      scale_color_gradient2(low="blue",high="red", mid="grey", 
                            limits=c(-maxAbsNES, maxAbsNES)) +
      theme_classic(base_size=12) + scale_size_continuous(limits=c(0,maxLogFDR)) + 
      ggplot2::labs(color = "NES", size = "-log(FDR)",title=m) +  
      theme(axis.title=element_blank(), plot.title=element_text(hjust=0.5, face="bold"),
            axis.text.x = element_text(angle=45, vjust=1, hjust=1)) +
      geom_point(data = subset(dot.df, sig), col = "black", stroke = 1.5, shape = 21)
    ggplot2::ggsave(paste0(fname,"_inhibitors.pdf"), dot.plot, width=(min.width-1), height=customHeight)
    ggplot2::ggsave(paste0(fname,"_inhibitors_wider.pdf"), dot.plot, width=min.width, height=customHeight)
    ggplot2::ggsave(paste0(fname,"_inhibitors_evenWider.pdf"), dot.plot, width=(min.width+1), height=customHeight)
  } else if (customHeight<49) {
    ggplot2::ggsave(paste0(fname,"_wider.pdf"), dot.plot, width=(min.width+1), height=customHeight)
    ggplot2::ggsave(paste0(fname,"_widerTaller.pdf"), dot.plot, width=(min.width+1), height=(customHeight+1))
  }
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

sens <- read.csv(synapser::synGet("syn65672258")$path)
sens$sig <- FALSE
sens[sens$p_value <= 0.05 & sens$FDR_q_value <= 0.25,]$sig <- TRUE
sens$Timepoint <- sens$Time
sens$Timepoint <- factor(sens$Timepoint, levels=c("8h","24h"))
sens[sens$FDR_q_value == 0,]$minusLogFDR <- 4
sens[sens$Drug_set == "Bcr.Abl kinase inhibitor",]$Drug_set <- "Bcr-Abl kinase inhibitor"

gsea <- read.csv(synapser::synGet("syn65928903")$path)
gsea$Timepoint <- gsea$Time
gsea$Timepoint <- factor(gsea$Timepoint, levels=c("8h", "24h"))
gsea$individualID <- gsea$MPNST
gsea$Drug <- gsea$DrugTreatment
gsea$signif <- gsea$sig
gsea$feature <- sub("KEGG_","",gsea$Feature_set)
gsea[gsea$FDR_q_value==0,]$minusLogFDR <- 4
keepCols2 <- c("Drug_set", "DrugTreatment","NES","minusLogFDR","sig","Ranking","Timepoint")

tf <- read.csv(synapser::synGet("syn64420968")$path)
tf$feature <- tf$source
tf$minusLogFDR <- -log10(tf$padj)

maxAbsNESVals <- list("Differential Expression" = max(abs(diffexp$log2FoldChange), na.rm=TRUE),
                      "Gene Set Enrichment" = max(abs(gsea$NES), na.rm=TRUE),
                      "Transcription Factor Enrichment" = max(abs(tf$score), na.rm=TRUE))
maxLogFDRVals <- list("Differential Expression" = max(diffexp$minusLogFDR, na.rm=TRUE),
                      "Gene Set Enrichment" = max(gsea$minusLogFDR, na.rm=TRUE),
                      "Transcription Factor Enrichment" = max(tf$minusLogFDR, na.rm=TRUE))

#### find drug targets ####
# for each drugTreatment MOA
drug.info <- list("Palbociclib" = "CDK inhibitor", "Ribociclib" = "CDK inhibitor",
                  "Trabectidin" = "Chemotherapy", "Ifosfamide" = "DNA alkylating agent",
                  "Decitabine" = "DNMT inhibitor", "Vorinostat" = "HDAC inhibitor",
                  "Mirdametinib" = "MEK inhibitor", "Selumetinib" = "MEK inhibitor",
                  "Trametinib" = "MEK inhibitor", "Capmatinib" = "MET inhibitor",
                  "Olaparib" = "PARP inhibitor", "RMC4630" = "SHP2 inhibitor",
                  "TNO155" = "SHP2 inhibitor", "Doxorubicin" = "TOP inhibitor",
                  "Irinotecan" = "TOP inhibitor", "Verteporfin" = "YAP inhibitor")
# also have version matching PRISM drug moa sets
drug.info2 <- list("Palbociclib" = "CDK inhibitor", "Ribociclib" = "CDK inhibitor",
                  "Trabectidin" = "Chemotherapy", "Ifosfamide" = "DNA alkylating agent",
                  "Decitabine" = "DNA methyltransferase inhibitor", "Vorinostat" = "HDAC inhibitor",
                  "Mirdametinib" = "MEK inhibitor", "Selumetinib" = "MEK inhibitor",
                  "Trametinib" = "MEK inhibitor", "Capmatinib" = "MET inhibitor",
                  "Olaparib" = "PARP inhibitor", "RMC4630" = "SHP2 inhibitor",
                  "TNO155" = "SHP2 inhibitor", "Doxorubicin" = "Topoisomerase inhibitor",
                  "Irinotecan" = "Topoisomerase inhibitor", "Verteporfin" = "YAP inhibitor")

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

maxAbsNES <- ceiling(max(abs(na.omit(diffexp[diffexp$hgnc_symbol %in% target.diffexp$hgnc_symbol &
                                               diffexp$Drug %in% target.diffexp$Drug,]$log2FoldChange)))) # 7
maxAbsP <- ceiling(max(-log10(na.omit(diffexp[diffexp$hgnc_symbol %in% target.diffexp$hgnc_symbol &
                                                diffexp$Drug %in% target.diffexp$Drug,]$padj)))) # 22

# plot of targets per MOA
incl.drugs <- unique(target.diffexp$Drug)
moa <- sort(unique(unlist(drug.info[incl.drugs]))) # 8
target.plots <- NULL
for (i in moa) {
  cat(i,"\n")
  temp.drugs <- names(drug.info[drug.info == i])
  temp.targets <- unique(target.diffexp[target.diffexp$Drug %in% temp.drugs,]$hgnc_symbol)
  temp.diffexp <- diffexp[diffexp$hgnc_symbol %in% temp.targets & diffexp$Drug %in% temp.drugs,]
  geneOrder <- unique(temp.diffexp[order(temp.diffexp$log2FoldChange, decreasing=TRUE),]$hgnc_symbol)
  tempMaxAbsNES <- ceiling(max(abs(na.omit(temp.diffexp$log2FoldChange))))
  tempMaxAbsP <- ceiling(max(-log10(na.omit(temp.diffexp$padj))))
  temp.plot <- ggplot2::ggplot(temp.diffexp, aes(x=Timepoint, y=hgnc_symbol, color=log2FoldChange, size=-log10(padj))) +
    facet_grid(individualID~Drug) + geom_point() + theme_classic(base_size=12) + 
    theme(axis.title=element_blank(),axis.text.x=element_text(angle=45, vjust=1, hjust=1), plot.title=element_text(hjust=0.5, face="bold")) + 
    labs(title=i, color="Log2FC", size="-Log(FDR)") + ggplot2::scale_y_discrete(limits = geneOrder) +
    scale_size_continuous(limits=c(0,tempMaxAbsP)) +
    scale_color_gradient2(low="blue",high="red", mid="grey", 
                          limits=c(-tempMaxAbsNES, tempMaxAbsNES)) +
    geom_point(data = subset(temp.diffexp, signif), col = "black", stroke = 1.5, shape = 21)
  tempWidth <- ifelse(length(temp.drugs)<2, 3, 2*length(temp.drugs)) # each moa has 1-3 drugs represented
  tempHeight <- ifelse(0.5*length(temp.targets)<2, 2.5, 0.5*length(temp.targets)) # each moa has 1-9 targets measured
  ggsave(paste0(i,"_target_dotPlot.pdf"), temp.plot, width=tempWidth, height=tempHeight)
  if (is.null(target.plots)) {
    target.plots <- temp.plot +
      scale_size_continuous(limits=c(0,maxAbsP)) +
      scale_color_gradient2(low="blue",high="red", mid="grey", 
                            limits=c(-maxAbsNES, maxAbsNES))
  } else {
    target.plots <- target.plots + 
      (temp.plot + scale_size_continuous(limits=c(0,maxAbsP)) + 
         scale_color_gradient2(low="blue",high="red", mid="grey", limits=c(-maxAbsNES, maxAbsNES))) + 
      plot_layout(guides="collect")
  }
}
target.plots <- target.plots + plot_layout(nrow=2)
ggsave("target_dotPlot.pdf", target.plots, width=9, height=8)
ggsave("target_dotPlot_wider.pdf", target.plots, width=12, height=8)
ggsave("target_dotPlot_widest.pdf", target.plots, width=16, height=8)

#### use results to predict synergy ####
# make plots for each moa
inputs <- list("Differential Expression" = diffexp, 
               "Gene Set Enrichment" = gsea,
               "Transcription Factor Enrichment" = tf)
maxAbsNES <- max(abs(sens$NES)) # 4.52
maxLogFDR <- max(sens$minusLogFDR) # 4
for (i in unique(unlist(drug.info2))) {
  # if samples are resistant to treatment, filter for sensitivity
  # else filter for resistance if they are still sensitive to the given treatment
  
  # likewise, if drug is similar to given MOA, filter for similarity
  # else filter for opposite behavior
  temp.drugs <- names(drug.info2[drug.info2 == i])
  mpnst <- na.omit(unique(sens$MPNST))
  timepoints <- levels(sens$Timepoint)
  keepCols <- c("MPNST", "Timepoint", "DrugTreatment", 
                "Drug_set", "NES", "minusLogFDR", "sig", "Ranking")
  for (m in mpnst) {
    venn.times <- list()
    venn.sens <- list()
    venn.plot.times <- NULL
    dmea.df <- data.frame()
    
    strict.venn.times <- list()
    strict.venn.sens <- list()
    strict.venn.plot.times <- NULL
    strict.dmea.df <- data.frame()
    for (t in timepoints) {
      venn.list <- list()
      strict.venn.list <- list()
      for (j in temp.drugs) {
        # get sig sensitivity
        temp.sens <- na.omit(sens[sens$MPNST == m & sens$Timepoint==t &
                                    sens$DrugTreatment == j,])
        given.sens <- mean(temp.sens[temp.sens$Drug_set == i & temp.sens$sig,]$NES, na.rm=TRUE)
        sens.dot.df <- na.omit(temp.sens[temp.sens$Drug_set %in% 
                                           unique(na.omit(temp.sens[temp.sens$sig,]$Drug_set)),])
        if (is.numeric(given.sens) & given.sens != "NaN") {
          strict.sens <- na.omit(sens[sens$MPNST == m & sens$Timepoint==t &
                                      sens$DrugTreatment == j & 
                                        sens$NES/given.sens < 0,])
          venn.color <- ifelse(given.sens < 0, "blue", "red")
        } else {
          strict.sens <- na.omit(sens[sens$MPNST == m & sens$Timepoint==t &
                                        sens$DrugTreatment == j & 
                                        sens$NES < 0,])
          venn.color <- "blue"
        }
        strict.sens.dot.df <- na.omit(temp.sens[temp.sens$Drug_set %in% 
                                                  unique(na.omit(strict.sens[strict.sens$sig,]$Drug_set)),])
        
        if (nrow(sens.dot.df) > 0) {
          if (any(sens.dot.df$minusLogFDR == "Inf")) {
            sens.dot.df[sens.dot.df$minusLogFDR == "Inf",]$minusLogFDR <- 4
          }
          sens.dot.df$Ranking <- "Sensitivity"
          dmea.df <- rbind(dmea.df, sens.dot.df[,keepCols])
          venn.list[[paste0("Sensitivity:\n",j)]] <- unique(na.omit(sens.dot.df[sens.dot.df$sig,]$Drug_set))
        }
        
        if (nrow(strict.sens.dot.df) > 0) {
          if (any(strict.sens.dot.df$minusLogFDR == "Inf")) {
            strict.sens.dot.df[strict.sens.dot.df$minusLogFDR == "Inf",]$minusLogFDR <- 4
          }
          strict.sens.dot.df$Ranking <- "Sensitivity"
          strict.dmea.df <- rbind(strict.dmea.df, strict.sens.dot.df[,keepCols])
          strict.venn.list[[paste0("Sensitivity:\n",j)]] <- unique(na.omit(strict.sens.dot.df[strict.sens.dot.df$sig,]$Drug_set))
        }
      }
        
      if (nrow(dmea.df) > 0) {
        # venn diagrams
        venn.list <- venn.list[sort(names(venn.list))]
        if (length(venn.list) > 1 & 
            length(venn.list) < 5) {
          venn.plot <- wrap_elements(
            ggvenn::ggvenn(venn.list, show_percentage = FALSE, 
                           set_name_size = 1, text_size = 4) + 
              plot_annotation(title=t, theme=theme(plot.title=element_text(hjust=0.5, face="bold"))))
        } else {venn.plot <- list()}
        venn.sens[[t]] <- venn.list[grepl("Sensitivity",names(venn.list))]
        venn.times[[t]] <- venn.list
        
        if (is.null(venn.plot.times) & length(venn.plot) > 0) {
          venn.plot.times <- venn.plot
        } else if (length(venn.plot) > 0) {
          venn.plot.times <- venn.plot.times / venn.plot
        }
      }
      
      if (nrow(strict.dmea.df) > 0) {
        # venn diagrams
        strict.venn.list <- strict.venn.list[sort(names(strict.venn.list))]
        if (length(strict.venn.list) > 1 & 
            length(strict.venn.list) < 5) {
          strict.venn.plot <- wrap_elements(
            ggvenn::ggvenn(strict.venn.list, show_percentage = FALSE, 
                           fill_color = rep(venn.color,length(strict.venn.list)),
                           set_name_size = 1, text_size = 4) + 
              plot_annotation(title=t, theme=theme(plot.title=element_text(hjust=0.5, face="bold"))))
        } else {strict.venn.plot <- list()}
        strict.venn.sens[[t]] <- strict.venn.list[grepl("Sensitivity",names(venn.list))]
        strict.venn.times[[t]] <- strict.venn.list
        
        if (is.null(strict.venn.plot.times) & length(strict.venn.plot) > 0) {
          strict.venn.plot.times <- strict.venn.plot
        } else if (length(strict.venn.plot) > 0) {
          strict.venn.plot.times <- strict.venn.plot.times / strict.venn.plot
        }
      }
    } 
    if (length(venn.list) > 0) {
      ggplot2::ggsave(paste0(m,"_",i,"_MOA_vennDiagram.pdf"), venn.plot.times, width=3, height=5) 
    }
    if (length(strict.venn.list) > 0) {
      ggplot2::ggsave(paste0(m,"_",i,"_oppNES_MOA_vennDiagram.pdf"), strict.venn.plot.times, width=3, height=5) 
    }
    
    # dot plots
    sigMOAs <- unique(unlist(venn.sens))
    if (length(sigMOAs) > 0) {
      mean.dmea.df <- plyr::ddply(dmea.df, .(Drug_set), summarize,
                                  absNES = mean(abs(NES), na.rm=TRUE))
      sigOrder <- mean.dmea.df[order(mean.dmea.df$absNES),]$Drug_set
      sigOrder <- sigOrder[sigOrder %in% sigMOAs]
      dot.df <- na.omit(dmea.df)
      DMEAdotPlots(dot.df, y.order = sigOrder, title=m, 
                   fname = paste0(m,"_",i,"_DMEA_dotPlot"), 
                   maxAbsNES = maxAbsNES, maxLogFDR = maxLogFDR)
      
      # also filter for opposite of given MOA sensitivity OR NES < 0 if NA
      if (is.numeric(given.sens) & given.sens !="NaN") {
        sig.dot.df <- dot.df[dot.df$NES/given.sens < 0 & dot.df$sig,]
      } else {
        sig.dot.df <- dot.df[dot.df$NES<0 & dot.df$sig,]
      }
      sigOrder <- sigOrder[sigOrder %in% sig.dot.df$Drug_set]
      dot.df <- dot.df[dot.df$Drug_set %in% sigOrder,]
      DMEAdotPlots(dot.df, y.order = sigOrder, title=m, 
                    fname = paste0(m,"_",i,"_oppNES_DMEA_dotPlot"), 
                    maxAbsNES = maxAbsNES, maxLogFDR = maxLogFDR)
      
      # also filter for shared opposite of given MOA sensitivity OR NES < 0 if NA
      strictSigShared <- c()
      for (t in names(strict.venn.sens)) {
        strictSigShared <- unique(c(strictSigShared, Reduce(intersect, strict.venn.sens[[t]])))
      }
      if (length(strictSigShared) > 0) {
        sigOrder <- sigOrder[sigOrder %in% strictSigShared]
        dot.df <- dot.df[dot.df$Drug_set %in% sigOrder,]
        DMEAdotPlots(dot.df, y.order = sigOrder, title=m, 
                     fname = paste0(m,"_",i,"_oppNESShared_DMEA_dotPlot"), 
                     maxAbsNES = maxAbsNES, maxLogFDR = maxLogFDR) 
      }
      
      # look at shared diffexp, GSEA, TF enrichment
      for (otherMOA in sigMOAs[sigMOAs != i]) {
        other.drugs <- names(drug.info2[drug.info2 == otherMOA])
        if (length(other.drugs) > 0) {
          both.drugs <- c(temp.drugs, other.drugs)
          for (temp.input in names(inputs)) {
            input <- inputs[[temp.input]]
            sigList <- list()
            for (drug in both.drugs) {
              if (nrow(input[input$Drug == drug & input$individualID == m & input$signif,]) > 0) {
                sigList[[drug]] <- na.omit(unique(input[input$Drug == drug & input$individualID == m & input$signif,]$feature)) 
              }
            }
            
            if (length(sigList) > 1) {
              if (length(sigList) > 4) {
                # convert to matrix
                myPlot = ComplexHeatmap::make_comb_mat(sigList)
                #ComplexHeatmap::UpSet(m, right_annotation = moaAnnotation)
                pdf(paste0(m,"_",i,"_",otherMOA,"_sig_",temp.input,"_upsetPlot.pdf"))
                up <- ComplexHeatmap::UpSet(myPlot)
                ComplexHeatmap::draw(up)
                dev.off() 
              } else if (length(sigList) > 1) {
                venn.plot <- wrap_elements(
                  ggvenn::ggvenn(sigList, show_percentage = FALSE, 
                                 set_name_size = 3, text_size = 4) + 
                    plot_annotation(title=m, theme=theme(plot.title=element_text(hjust=0.5, face="bold"))))
                ggplot2::ggsave(paste0(m,"_",i,"_",otherMOA, "_sig_",temp.input,"_vennDiagram.pdf"), venn.plot, width=3, height=5) 
              }
              
              # also make dot plot
              overlapSig <- unlist(sigList[temp.drugs])[unlist(sigList[temp.drugs]) %in% unlist(sigList[other.drugs])]
              if (length(overlapSig) > 0) {
                dot.df <- input[input$Drug %in% both.drugs & input$individualID == m & input$feature %in% overlapSig,]
                
                # change colnames so we don't have to make another custom dotPlot function
                dot.df$DrugTreatment <- dot.df$Drug
                dot.df$Drug_set <- dot.df$feature
                if (temp.input == "Differential Expression") {
                  dot.df$NES <- dot.df$log2FoldChange
                  dot.df$minusLogFDR <- -log10(dot.df$padj)
                } else if (temp.input == "Transcription Factor Enrichment") {
                  dot.df$NES <- dot.df$score
                  dot.df$minusLogFDR <- -log10(dot.df$padj)
                }
                dot.df$Ranking <- temp.input
                dot.df$sig <- dot.df$signif
                mean.df <- plyr::ddply(dot.df, .(Drug_set), summarize,
                                       absNES = mean(abs(NES), na.rm=TRUE))
                sigOrder <- mean.df[order(mean.df$absNES),]$Drug_set
                if (length(sigOrder)/5 < 49) {
                  if (temp.input == "Gene Set Enrichment") {
                    min.width=7
                  } else { min.width = 5 }
                  DMEAdotPlots(dot.df, y.order = sigOrder, title=m, min.width = min.width,
                               fname = paste0(m,"_",i,"_",otherMOA,"_",temp.input,"_dotPlot"), 
                               maxAbsNES = maxAbsNESVals[[temp.input]], maxLogFDR = maxLogFDRVals[[temp.input]]) 
                } 
                
                # also run correlation(s)
                dot.df$NES <- as.numeric(dot.df$NES)
                corr.input <- na.omit(dot.df[,c("feature","Drug","NES")])
                Log2FC.df <- reshape2::dcast(corr.input, feature ~ Drug, fun.aggregate=mean,
                                             value.var = "NES")
                
                # convert from data frame to numeric matrix
                rownames(Log2FC.df) <- Log2FC.df$feature
                Log2FC.mat <- as.matrix(Log2FC.df[, 2:ncol(Log2FC.df)])
                
                # create correlation matrix
                corr.mat <- stats::cor(Log2FC.mat)
                write.csv(corr.mat,paste0(m,"_",i,"_",otherMOA,"_sig_",temp.input,"_PearsonCorr.csv"))
                
                # plot correlation matrix
                corr.mat.plot <- ggcorrplot::ggcorrplot(corr.mat)
                ggsave(paste0(m,"_",i,"_",otherMOA,"_sig_",temp.input,"_PearsonCorr.pdf"), corr.mat.plot, width=4, height=4)
                
                # check for known targets
                temp.targets <- unique(target.diffexp[target.diffexp$Drug %in% temp.drugs,]$hgnc_symbol)
                if (any(temp.targets %in% overlapSig)) {
                  dot.df <- input[input$Drug %in% both.drugs & input$feature %in% temp.targets & input$individualID == m,]
                  # change colnames so we don't have to make another custom dotPlot function
                  dot.df$DrugTreatment <- dot.df$Drug
                  dot.df$Drug_set <- dot.df$feature
                  if (temp.input == "Differential Expression") {
                    dot.df$NES <- dot.df$log2FoldChange
                    dot.df$minusLogFDR <- -log10(dot.df$padj)
                  } else if (temp.input == "Transcription Factor Enrichment") {
                    dot.df$NES <- dot.df$score
                    dot.df$minusLogFDR <- -log10(dot.df$padj)
                  }
                  dot.df$Ranking <- temp.input
                  dot.df$sig <- dot.df$signif
                  mean.df <- plyr::ddply(dot.df, .(Drug_set), summarize,
                                         absNES = mean(abs(NES), na.rm=TRUE))
                  sigOrder <- mean.df[order(mean.df$absNES),]$Drug_set
                  if (length(sigOrder)/5 < 49) {
                    if (temp.input == "Gene Set Enrichment") {
                      min.width=7
                    } else { min.width = 5 }
                    DMEAdotPlots(dot.df, y.order = sigOrder, title=m, min.width = min.width,
                                 fname = paste0(m,"_",i,"_",otherMOA,"_",temp.input,"_target_dotPlot"), 
                                 maxAbsNES = maxAbsNESVals[[temp.input]], maxLogFDR = maxLogFDRVals[[temp.input]]) 
                  }
                }
              }  
            }
          }
        }
      }
    }
  }
}

# plot top DEGs
dir.create("diffexp_dotPlots")
setwd("diffexp_dotPlots")
maxAbsNES <- max(abs(diffexp$log2FoldChange), na.rm=TRUE)
maxLogFDR <- max(-log10(diffexp$padj), na.rm=TRUE)
for (m in mpnst) {
  for (t in timepoints) {
    for (d in tested.drugs) {
      dot.df <- diffexp[diffexp$individualID == m & diffexp$Timepoint == t & diffexp$Drug == d,]
      sig.dot.df <- na.omit(dot.df[dot.df$signif,] %>% slice_max(abs(log2FoldChange),n=10))
      geneOrder <- unique(sig.dot.df[order(sig.dot.df$log2FoldChange),]$hgnc_symbol)
      dot.df <- na.omit(dot.df[dot.df$hgnc_symbol %in% sig.dot.df$hgnc_symbol,])
      if (nrow(sig.dot.df) > 0) {
        dot.df$DrugTreatment <- dot.df$Drug
        dot.df$Drug_set <- dot.df$hgnc_symbol
        dot.df$NES <- dot.df$log2FoldChange
        dot.df$minusLogFDR <- -log10(dot.df$padj)
        dot.df$Ranking <- "Differential Expression"
        dot.df$sig <- dot.df$signif
        DMEAdotPlots(dot.df, y.order=geneOrder, title=m,
                     fname=paste0(m,"_",d,"_",t,"_top10diffexp_dotPlot"), 
                     maxAbsNES=maxAbsNES, maxLogFDR = maxLogFDR) 
      } else {
        warning("no DEGs for",m,d,t)
      }
    }
  }
}

#### TF network centrality ####
# what are the most central proteins in each MOA's TF network(s)?
centrality <- data.frame()
mpnst <- c("JH-2-002","MN-2")
timepoints <- c("8h","24h")
tf.path <- "/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/MPNST-PDX-MT/RNAseq/synapse_based_code/TF_networks"
for (i in unique(unlist(drug.info2))) {
  temp.drugs <- names(drug.info2[drug.info2 == i])
  for (m in mpnst) {
    for (t in timepoints) {
      for (j in temp.drugs) {
        fname <- paste0(m,"_",j,"_",t)
        pos.fname <- file.path(tf.path,paste0(fname, "_positive_TF_centrality.csv"))
        neg.fname <- file.path(tf.path,paste0(fname, "_negative_TF_centrality.csv"))
        if (file.exists(pos.fname)) {
          pos.centrality <- read.csv(pos.fname)
          pos.centrality$Direction <- "positive"
          pos.centrality$Drug <- j
          pos.centrality$Timepoint <- t
          pos.centrality$individualID <- m
          pos.centrality$moa <- i
          centrality <- rbind(centrality, pos.centrality)
        }
        if (file.exists(neg.fname)) {
          neg.centrality <- read.csv(neg.fname) 
          neg.centrality$Direction <- "Negative"
          neg.centrality$Drug <- j
          neg.centrality$Timepoint <- t
          neg.centrality$individualID <- m
          neg.centrality$moa <- i
          centrality <- rbind(centrality, neg.centrality)
        }
      }
    }
  }
}
write.csv(centrality, "TF_centrality.csv", row.names=FALSE)

mean.centrality <- plyr::ddply(centrality,.(individualID, moa, name), summarize,
                               meanCentrality=mean(eigen_centrality, na.rm=TRUE),
                               drugs = paste0(sort(unique(Drug)), collapse=", "),
                               directions=paste0(sort(unique(Direction)), collapse=", "),
                               times = paste0(sort(unique(Timepoint)), collapse=", "))
write.csv(mean.centrality, "TF_meanCentrality.csv", row.names=FALSE)

# #### look at correlations ####
corr <- read.csv("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/MPNST-PDX-MT/RNAseq/synapse_based_code/degCorr.csv")
rownames(corr) <- corr$X
corr$X <- NULL
setwd("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/MPNST-PDX-MT/RNAseq/synapse_based_code")
dir.create("diffexpCorrelations")
setwd("diffexpCorrelations")
library(plyr);library(dplyr);library(tidyr)
for (i in unique(unlist(drug.info2))) {
  temp.drugs <- names(drug.info2[drug.info2 == i])
  moa.cols <- c()
  for (d in temp.drugs) {
    moa.cols <- c(colnames(corr)[grepl(d,colnames(corr))],moa.cols) 
  }
  moa.corr <- corr[,moa.cols]
  moa.corr$otherDrug <- rownames(moa.corr)
  moa.corr.df <- reshape2::melt(moa.corr)
  moa.corr.df <- moa.corr.df %>% tidyr::separate_wider_delim(variable,"_",names=c("MPNST","Timepoint","Drug"))
  moa.corr.df <- moa.corr.df %>% tidyr::separate_wider_delim(otherDrug,"_",names=c("MPNST2","Timepoint2","Drug2"))
  moa.corr.df$MPNST <- gsub("[.]","-",moa.corr.df$MPNST)
  moa.corr.df <- moa.corr.df[moa.corr.df$MPNST == moa.corr.df$MPNST2 & moa.corr.df$Timepoint==moa.corr.df$Timepoint2,]
  # mean.moa.corr.df <- plyr::ddply(moa.corr.df,.(MPNST, Drug, Drug2), summarize,
  #                                 meanPearsonEst = mean(value, na.rm=TRUE),
  #                                 medianPearsonEst = median(value,na.rm=TRUE),
  #                                 sdPearsonEst = sd(value, na.rm=TRUE))
  moa.corr.df$Timepoint <- factor(moa.corr.df$Timepoint, levels=c("8h","24h"))
  moa.corr.plot <- ggplot(moa.corr.df, aes(x=Drug2, y=value)) + geom_violin(alpha=0) +
    geom_point(aes(color=Timepoint, shape=Drug)) +
    geom_boxplot(width=0.2, alpha = 0) +facet_wrap(.~MPNST)+theme_classic()+
    scale_x_discrete(limits=names(drug.info2)) + ylab("Pearson Correlation") +
    theme(axis.title.x=element_blank(),axis.text.x=element_text(angle=45, hjust=1, vjust=1),
          plot.title=element_text(face="bold", hjust=0.5)) +
    ggtitle(i) + scale_y_continuous(limits=c(0,1))
  ggsave(paste0(i,"_diffexp_correlations_v2.pdf"),moa.corr.plot,width=7,height=2)
  moa.corr.plot <- ggplot(moa.corr.df, aes(x=Drug2, y=value)) + geom_violin(alpha=0) +
    geom_point(aes(color=Timepoint, shape=Drug)) +
    geom_boxplot(width=0.2, alpha = 0) +facet_wrap(.~MPNST)+theme_classic()+
    scale_x_discrete(limits=names(drug.info2)) + ylab("Pearson Correlation") +
    theme(axis.title.x=element_blank(),axis.text.x=element_text(angle=45, hjust=1, vjust=1),
          plot.title=element_text(face="bold", hjust=0.5)) +
    ggtitle(i) + scale_y_continuous(limits=c(-1,1))
  ggsave(paste0(i,"_diffexp_correlations_v2.pdf"),moa.corr.plot,width=7,height=4)
}

