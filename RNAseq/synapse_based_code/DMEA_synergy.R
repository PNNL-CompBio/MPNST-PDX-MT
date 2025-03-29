# compile DMEA for potential synergy
library(plyr);library(dplyr);library(synapser);library(data.table);library(ggplot2);library(patchwork)
synapser::synLogin()
setwd("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/MPNST-PDX-MT/RNAseq/synapse_based_code")
dir.create("DMEA_synergy")
setwd("DMEA_synergy")

# load results
diffexp <- read.csv(synapser::synGet('syn64397743')$path)
diffexp$log2FoldChange <- as.numeric(diffexp$log2FoldChange)
diffexp$padj <- as.numeric(diffexp$padj)
diffexp$Timepoint <- factor(diffexp$Timepoint, levels=c("8h", "24h"))
sens <- read.csv(synapser::synGet("syn65672258")$path)
sens$sig <- FALSE
sens[sens$p_value <= 0.05 & sens$FDR_q_value <= 0.25,]$sig <- TRUE
colnames(sens)[1] <- "Timepoint" # change 'Time' to 'Timepoint'
sens$Timepoint <- factor(sens$Timepoint, levels=c("8h","24h"))
sim <- read.csv(synapser::synGet("syn65674473")$path)
sim$sig <- FALSE
sim[sim$p_value <= 0.05 & sim$FDR_q_value <= 0.25,]$sig <- TRUE
sim$Timepoint <- paste0(sim$Time, "h")
sim$Timepoint <- factor(sim$Timepoint, levels=c("8h","24h"))

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
                  "TNO155" = "SHP2 inhibitor", "Doxorubicin" = "TOP inhibitor",
                  "Irinotecan" = "topoisomerase inhibitor", "Verteporfin" = "YAP inhibitor")
#drug.target.info <- read.csv("https://github.com/BelindaBGarana/DMEA/raw/refs/heads/shiny-app/Inputs/PRISM_secondary-screen-replicate-treatment-info.csv")
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

#### find optimal timepoint for each drug ####
measured.targets <- list()
Drug <- names(drug.info)
time.info <- data.frame()
optTimes <- data.frame(Drug, Timepoint = NA)
for (i in names(drug.info)){
  temp.targets <- data.frame()
  temp.diffexp <- diffexp[diffexp$hgnc_symbol != "" & diffexp$Drug == i,]
  for (j in drug.targets[[i]]) {
    temp.targets <- rbind(temp.targets, na.omit(temp.diffexp[startsWith(j,temp.diffexp$hgnc_symbol),])) 
  }
  if (nrow(temp.targets) > 0) {
    measured.targets[[i]] <- temp.targets 
    mean.effect <- plyr::ddply(temp.targets, .(Drug, Timepoint), summarize,
                               meanLog2FC = mean(log2FoldChange, na.rm = TRUE))
    time.info <- rbind(time.info, mean.effect)
    optTimes[optTimes$Drug == i,]$Timepoint <- levels(mean.effect$Timepoint)[mean.effect[which.max(abs(mean.effect$meanLog2FC)),]$Timepoint]
  }
}
target.diffexp <- data.table::rbindlist(measured.targets, use.names=FALSE)
write.csv(target.diffexp, "target_diffexp.csv", row.names = FALSE)
write.csv(time.info, "meanTarget_diffexp.csv", row.names=FALSE)
write.csv(optTimes, "optimalTreatmentTimes.csv", row.names=FALSE)
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

#### use results from optimal times to predict synergy ####
# if no optimal time determined due to lack of known/measured targets,
# use timepoint with most diffexp genes
# else if no degs, then use 24h
opt.sens <- data.frame()
opt.sim <- data.frame()
for (i in names(drug.info)) {
  # find optimal time
  temp.time <- optTimes[optTimes$Drug == i,]$Timepoint
  drug.diffexp <- na.omit(diffexp[diffexp$Drug==i & diffexp$signif,])
  if (is.na(temp.time) & nrow(drug.diffexp) == 0) {
    temp.time <- "24h"
  } else if (is.na(temp.time)) {
    n.degs <- list()
    for (j in unique(drug.diffexp$Timepoint)) {
      n.degs[[j]] <- nrow(drug.diffexp[drug.diffexp$Timepoint == j,])
    }
    temp.time <- names(n.degs)[which.max(n.degs)]
  }
  
  # get results at optimal time
  temp.sens <- na.omit(sens[sens$Timepoint==temp.time & sens$DrugTreatment==i,])
  opt.sens <- rbind(opt.sens, temp.sens)
  temp.sim <- na.omit(sim[sim$Timepoint==temp.time & sim$DrugTreatment==i,])
  opt.sim <- rbind(opt.sim, temp.sim)
}

# make venn diagrams for each moa
for (i in unique(unlist(drug.info2))) {
  # if samples are resistant to treatment, filter for sensitivity
  # else filter for resistance if they are still sensitive to the given treatment
  
  # likewise, if drug is similar to given MOA, filter for similarity
  # else filter for opposite behavior
  
  venn.list <- list()
  strict.venn.list <- list()
  venn.colors <- list()
  temp.drugs <- names(drug.info2[drug.info2 == i])
  
  given.sens <- opt.sens[opt.sens$DrugTreatment %in% temp.drugs & 
                           opt.sens$Drug_set == i,]
  if (nrow(given.sens) > 0) {
    dir.given.sens <- mean(given.sens[given.sens$sig,]$NES, na.rm=TRUE)
    if (is.na(dir.given.sens)) {
      dir.given.sens <-  mean(given.sens$NES, na.rm=TRUE)
    } 
  } else {
    dir.given.sens <- NULL
  }
  
  given.sim <- opt.sim[opt.sim$DrugTreatment %in% temp.drugs,]
  if (nrow(given.sim) > 0) {
    dir.given.sim <- mean(given.sim[given.sim$sig,]$NES, na.rm=TRUE)
    if (is.na(dir.given.sim)) {
      dir.given.sim <-  mean(given.sim$NES, na.rm=TRUE)
    }
  } else {
    dir.given.sim <- NULL
  }
  
  # we expect dir.given.sens and dir.given.sim to be opposite
  if (dir.given.sens/dir.given.sim > 0) {
    warning("sensitivity and similarity do not match for ",i)
  }
  
  for (j in temp.drugs) {
    drug.sens <- na.omit(opt.sens[opt.sens$DrugTreatment == j,])
    drug.sim <- na.omit(opt.sim[opt.sim$DrugTreatment == j,])
    
    venn.list[[paste0("Sensitivity:\n",j)]] <- unique(na.omit(drug.sens[drug.sens$sig,]$Drug_set))
    venn.list[[paste0("Similarity:\n",j)]] <- unique(na.omit(drug.sim[drug.sim$sig,]$Drug_set))
    
    strict.venn.list[[paste0("Sensitivity:\n",j)]] <- unique(na.omit(drug.sens[drug.sens$sig & 
                                                                                 sign(drug.sens$NES)/sign(dir.given.sens) < 0,]$Drug_set))
    strict.venn.list[[paste0("Similarity:\n",j)]] <- unique(na.omit(drug.sim[drug.sim$sig & 
                                                                               sign(drug.sim$NES)/sign(dir.given.sim) > 0,]$Drug_set))
    
    #venn.colors[[paste0("Sensitivity:\n",j)]] <- elseif(sign(dir.given.sens) < 0, )
  }
  venn.list <- venn.list[sort(names(venn.list))]
  strcit.venn.list <- strict.venn.list[sort(names(strict.venn.list))]
  
  ggvenn::ggvenn(venn.list, show_percentage = FALSE, set_name_size = 1, text_size = 4, show_elements=FALSE)
  ggplot2::ggsave("MOA_sensitivityANDsimilarity_vennDiagram.pdf", width=3, height=5)
  ggvenn::ggvenn(venn.list, show_percentage = FALSE, set_name_size = 2.5, text_size = 2.5, show_elements=TRUE)
  ggplot2::ggsave("MOA_sensitivityANDsimilarity_vennDiagram_wElements.pdf", width=5, height=5)
  
  ggvenn::ggvenn(strict.venn.list, show_percentage = FALSE, set_name_size = 1, text_size = 4, show_elements=FALSE,
                 fill_color = venn.colors)
  ggplot2::ggsave("MOA_sensitivityANDsimilarity_vennDiagram_oppSensNES_sameSimNES.pdf", width=3, height=5)
  ggvenn::ggvenn(strict.venn.list, show_percentage = FALSE, set_name_size = 2.5, text_size = 2.5, show_elements=TRUE,
                 fill.color=venn.colors)
  ggplot2::ggsave("MOA_sensitivityANDsimilarity_vennDiagram_oppSensNES_sameSimNES_wElements.pdf", width=5, height=5)
  # adjust venn colors bc not all moas have 2 drugs
  
}
venn.list <- list("Sensitivity:\nZhang et al" = sens[sens$type=="GSE183307" & sens$sig,]$Drug_set,
                  "Sensitivity:\nWilliams et al" = sens[sens$type=="GSE262030" & sens$sig,]$Drug_set,
                  "Similarity:\nWilliams et al" = sim[sim$Source=="GSE262030" & sim$sig,]$Drug_set,
                  "Similarity:\nZhang et al" = sim[sim$Source=="GSE183307" & sim$sig,]$Drug_set)
ggvenn::ggvenn(venn.list, show_percentage = FALSE, set_name_size = 1, text_size = 4, show_elements=FALSE)
ggplot2::ggsave("MOA_sensitivityANDsimilarity_vennDiagram.pdf", width=3, height=5)
ggvenn::ggvenn(venn.list, show_percentage = FALSE, set_name_size = 2.5, text_size = 2.5, show_elements=TRUE)
ggplot2::ggsave("MOA_sensitivityANDsimilarity_vennDiagram_wElements.pdf", width=5, height=5)

venn.list <- list("Sensitivity:\nGSE183307" = sens[sens$type=="GSE183307" & sens$sig & sens$NES<0,]$Drug_set,
                  "Sensitivity:\nGSE262030" = sens[sens$type=="GSE262030" & sens$sig & sens$NES<0,]$Drug_set,
                  "Similarity:\nGSE183307" = sim[sim$Source=="GSE183307" & sim$sig,]$Drug_set,
                  "Similarity:\nGSE262030" = sim[sim$Source=="GSE262030" & sim$sig,]$Drug_set)
ggvenn::ggvenn(venn.list, show_percentage = FALSE, set_name_size = 2.5, text_size = 5, show_elements=FALSE)
ggplot2::ggsave("MOA_sensitivityANDsimilarity_vennDiagram_negSensNES.pdf", width=5, height=5)
ggvenn::ggvenn(venn.list, show_percentage = FALSE, set_name_size = 2.5, text_size = 2.5, show_elements=TRUE)
ggplot2::ggsave("MOA_sensitivityANDsimilarity_vennDiagram_negSensNES_wElements.pdf", width=5, height=5)

venn.list <- list("Sensitivity:\nGSE183307" = sens[sens$type=="GSE183307" & sens$sig & sens$NES<0,]$Drug_set,
                  "Sensitivity:\nGSE262030" = sens[sens$type=="GSE262030" & sens$sig & sens$NES<0,]$Drug_set,
                  "Similarity:\nGSE262030" = sim[sim$Source=="GSE262030" & sim$sig & sim$NES>0,]$Drug_set,
                  "Similarity:\nGSE183307" = sim[sim$Source=="GSE183307" & sim$sig & sim$NES>0,]$Drug_set)
ggvenn::ggvenn(venn.list, show_percentage = FALSE, set_name_size = 1, text_size = 4, show_elements=FALSE)
ggplot2::ggsave("MOA_sensitivityANDsimilarity_vennDiagram_negSensNES_posSimNES.pdf", width=3, height=5)
ggvenn::ggvenn(venn.list, show_percentage = FALSE, set_name_size = 2.5, text_size = 2.5, show_elements=TRUE)
ggplot2::ggsave("MOA_sensitivityANDsimilarity_vennDiagram_negSensNES_posSimNES_wElements.pdf", width=5, height=5)

# create combined dot plot of DMEA results
colnames(sens)[1] <- "Source"
sens[sens$minusLogFDR == "Inf",]$minusLogFDR <- 4
sens$Ranking <- "Sensitivity"
sim$minusLogFDR <- 4
sim[sim$FDR_q_value !=0,]$minusLogFDR <- -log10(sim[sim$FDR_q_value !=0,]$FDR_q_value)
sim$Ranking <- "Similarity"
keepCols <- c("Source", "Drug_set", "NES", "minusLogFDR", "sig", "Ranking")
dmea.df <- rbind(sens[,keepCols], sim[,keepCols])
sigMOAs <- unique(c(sens[sens$sig,]$Drug_set, sim[sim$sig,]$Drug_set)) # 46
sharedMOAs <- sens$Drug_set[sens$Drug_set %in% sim$Drug_set] # 52
sigSharedMOAs <- sigMOAs[sigMOAs %in% sharedMOAs] # 15
mean.dmea.df <- plyr::ddply(dmea.df, .(Drug_set), summarize,
                            absNES = mean(abs(NES), na.rm=TRUE))
sigOrder <- mean.dmea.df[order(mean.dmea.df$absNES),]$Drug_set
sigOrder <- sigOrder[sigOrder %in% sigSharedMOAs]
dot.df <- dmea.df[dmea.df$Drug_set %in% sigSharedMOAs,]
library(ggplot2)
maxAbsNES <- max(abs(dot.df$NES))
dot.plot <- ggplot2::ggplot(
  dot.df,
  ggplot2::aes(
    x = Source, y = Drug_set, color = NES,
    size = minusLogFDR
  )
) + facet_wrap(. ~ Ranking) +
  ggplot2::geom_point() +
  ggplot2::scale_y_discrete(limits = sigOrder) +
  scale_color_gradient2(low="blue",high="red", mid="grey", 
                        limits=c(-maxAbsNES, maxAbsNES)) +
  theme_classic(base_size=12) + 
  theme(axis.title=element_blank(), 
        axis.text.x = element_text(angle=45, vjust=1, hjust=1)) +
  ggplot2::labs(color = "NES", size = "-log(FDR)") + 
  geom_point(data = subset(dot.df, sig), col = "black", stroke = 1.5, shape = 21)
dot.plot
ggplot2::ggsave("Compiled_DMEA_dotPlot_sensitivityANDsimilarity.pdf", width=4.2, height=4)

dot.df$Inhibitor <- sub(" inhibitor", "", dot.df$Drug_set)
sigOrder2 <- sub(" inhibitor", "", sigOrder)
dot.df[dot.df$Inhibitor == "Aurora kinase",]$Inhibitor <- "AURK"
sigOrder2 <- sub("Aurora kinase", "AURK", sigOrder2)
dot.plot <- ggplot2::ggplot(
  dot.df,
  ggplot2::aes(
    x = Source, y = Inhibitor, color = NES,
    size = minusLogFDR
  )
) + facet_wrap(. ~ Ranking) +
  ggplot2::geom_point() +
  ggplot2::scale_y_discrete(limits = sigOrder2) +
  scale_color_gradient2(low="blue",high="red", mid="grey", 
                        limits=c(-maxAbsNES, maxAbsNES)) +
  theme_classic(base_size=12) + 
  theme(axis.title=element_blank(), 
        axis.text.x = element_text(angle=45, vjust=1, hjust=1)) +
  ggplot2::labs(color = "NES", size = "-log(FDR)") + 
  geom_point(data = subset(dot.df, sig), col = "black", stroke = 1.5, shape = 21)
dot.plot
ggplot2::ggsave("Compiled_DMEA_dotPlot_sensitivityANDsimilarity_v2.pdf", width=3.3, height=4)

# just plot those neg in sens and pos in sim
sensMOAs <- dot.df[dot.df$NES < 0 & dot.df$Ranking == "Sensitivity",]$Inhibitor
simMOAs <- dot.df[dot.df$NES > 0 & dot.df$Ranking == "Similarity",]$Inhibitor
synMOAs <- simMOAs[simMOAs %in% sensMOAs] # 13
dot.df2 <- dot.df[dot.df$Inhibitor %in% synMOAs,]
sigOrder3 <- sigOrder2[sigOrder2 %in% synMOAs] # 7
dot.plot <- ggplot2::ggplot(
  dot.df2,
  ggplot2::aes(
    x = Source, y = Inhibitor, color = NES,
    size = minusLogFDR
  )
) + facet_wrap(. ~ Ranking) +
  ggplot2::geom_point() +
  ggplot2::scale_y_discrete(limits = sigOrder3) +
  scale_color_gradient2(low="blue",high="red", mid="grey", 
                        limits=c(-maxAbsNES, maxAbsNES)) +
  theme_classic(base_size=12) + 
  theme(axis.title=element_blank(), 
        axis.text.x = element_text(angle=45, vjust=1, hjust=1)) +
  ggplot2::labs(color = "NES", size = "-log(FDR)") + 
  geom_point(data = subset(dot.df2, sig), col = "black", stroke = 1.5, shape = 21)
dot.plot
ggplot2::ggsave("Compiled_DMEA_dotPlot_sensitivityANDsimilarity_v3.pdf", width=3.3, height=2.75)

# just plot those neg in sens and pos in sim - stricter
sensMOAs <- dot.df[dot.df$NES < 0 & dot.df$Ranking == "Sensitivity" & dot.df$sig,]$Inhibitor
simMOAs <- dot.df[dot.df$NES > 0 & dot.df$Ranking == "Similarity" & dot.df$sig,]$Inhibitor
synMOAs <- simMOAs[simMOAs %in% sensMOAs]
dot.df3 <- dot.df[dot.df$Inhibitor %in% synMOAs,]
sigOrder4 <- sigOrder2[sigOrder2 %in% synMOAs] # 3
dot.plot <- ggplot2::ggplot(
  dot.df3,
  ggplot2::aes(
    x = Source, y = Inhibitor, color = NES,
    size = minusLogFDR
  )
) + facet_wrap(. ~ Ranking) +
  ggplot2::geom_point() +
  ggplot2::scale_y_discrete(limits = sigOrder4) +
  scale_color_gradient2(low="blue",high="red", mid="grey", 
                        limits=c(-maxAbsNES, maxAbsNES)) +
  theme_classic(base_size=12) + 
  theme(axis.title=element_blank(), 
        axis.text.x = element_text(angle=45, vjust=1, hjust=1)) +
  ggplot2::labs(color = "NES", size = "-log(FDR)") + 
  geom_point(data = subset(dot.df3, sig), col = "black", stroke = 1.5, shape = 21) + 
  theme(legend.box = "horizontal")
# +
#theme(legend.position = "bottom"#, legend.box = "vertical"
#)
dot.plot
ggplot2::ggsave("Compiled_DMEA_dotPlot_sensitivityANDsimilarity_v4.pdf", width=4, height=2)
