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

# how many MOAs are evaluated for both sens and sim?
allSharedMOAs <- na.omit(unique(sens$Drug_set[sens$Drug_set %in% sim$Drug_set])) # 26

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
testedMOAsEvaluatable <- unique(unlist(drug.info2)[unlist(drug.info2) %in% allSharedMOAs])
testedMOAsNotEvaluatable <- unique(unlist(drug.info2)[!(unlist(drug.info2) %in% allSharedMOAs)]) # chemo, DNMTi, SHP2i, TOPi, YAPi
# in sim, decitabine is a "DNA inhibitor". Let's assume DNA inhibitor means DNMT inhibitor
#sim[sim$Drug_set == "DNA inhibitor",]$Drug_set <- "DNA methyltransferase inhibitor" # sens has "DNA methyltransferase inhibitor", "DNA inhibitor", "DNA polymerase inhibitor", "DNA synthesis inhibitor"

# other MOAs which are just named differently in sim vs. sens?
allSharedMOAs <- na.omit(unique(sens$Drug_set[tolower(sens$Drug_set) %in% tolower(sim$Drug_set)])) # 63! Much more than 26
MOAsUniqueToSens <- na.omit(unique(sens$Drug_set[!(tolower(sens$Drug_set) %in% tolower(sim$Drug_set))])) # 21
MOAsUniqueToSim <- na.omit(unique(sim$Drug_set[!(tolower(sim$Drug_set) %in% tolower(sens$Drug_set))])) # 128
sens[sens$Drug_set == "Bcr.Abl kinase inhibitor",]$Drug_set <- "Bcr-Abl inhibitor" # sens also has "Abl kinase inhibitor"?
sens[sens$Drug_set == "ALK tyrosine kinase receptor inhibitor",]$Drug_set <- "ALK inhibitor"
sens[sens$Drug_set == "glycogen synthase kinase inhibitor",]$Drug_set <- "GSK inhibitor"
sens[sens$Drug_set == "tubulin polymerization inhibitor",]$Drug_set <- "Tubulin inhibitor" # sens also has microtubule inhibitor
sens[sens$Drug_set == "topoisomerase inhibitor",]$Drug_set <- "Topoisomerase inhibitor"
sens[sens$Drug_set == "NFkB pathway inhibitor",]$Drug_set <- "NFKB inhibitor"
sens[sens$Drug_set == "PDGFR tyrosine kinase receptor inhibitor",]$Drug_set <- "PDGFR inhibitor"
sens[sens$Drug_set == "phosphatase inhibitor",]$Drug_set <- "Tyrosine phosphatase inhibitor"
sens[sens$Drug_set == "potassium channel blocker",]$Drug_set <- "Potassium channel antagonist"
sens[sens$Drug_set == "ribonucleotide reductase inhibitor",]$Drug_set <- "Ribonucleotide reductase inhibitor"
sim[sim$Drug_set == "Ribonucleoside reductase inhibitor",]$Drug_set <- "Ribonucleotide reductase inhibitor" # I think ribonucleotide makes more sense?
sens[sens$Drug_set == "sodium channel blocker",]$Drug_set <- "Sodium channel inhibitor"
sens[sens$Drug_set == "TGF beta receptor inhibitor",]$Drug_set <- "TGF-beta receptor inhibitor"
# sens has growth factor receptor inhibitor but sim has only EGFR/FGFR/PDGFR/VEGFR inhibitors (more specific)
# KIT inhibitor, RET tyrosine kinase inhibitor only in sens

allSharedMOAs <- sort(na.omit(unique(sens$Drug_set[tolower(sens$Drug_set) %in% tolower(sim$Drug_set)]))) # 74
allSharedMOAsSim <- sort(na.omit(unique(sim$Drug_set[tolower(sim$Drug_set) %in% tolower(sens$Drug_set)])))
MOAsUniqueToSens <- na.omit(unique(sens$Drug_set[!(tolower(sens$Drug_set) %in% tolower(sim$Drug_set))])) # 9
MOAsUniqueToSim <- na.omit(unique(sim$Drug_set[!(tolower(sim$Drug_set) %in% tolower(sens$Drug_set))])) # 117
moaMap <- data.frame(allSharedMOAs, allSharedMOAsSim)
for (i in 1:nrow(moaMap)) {
  sens[sens$Drug_set == moaMap$allSharedMOAs[i],]$Drug_set <- moaMap$allSharedMOAsSim[i]
}

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
#time.info <- data.frame()
optTimes <- data.frame(Drug, Timepoint = NA)
max.effect <- data.frame(Drug, Gene = NA, Log2FC = NA, p=NA, q=NA, Timepoint = NA, MPNST = NA)
for (i in names(drug.info)){
  temp.targets <- data.frame()
  temp.diffexp <- diffexp[diffexp$hgnc_symbol != "" & diffexp$Drug == i,]
  for (j in drug.targets[[i]]) {
    temp.targets <- rbind(temp.targets, na.omit(temp.diffexp[startsWith(j,temp.diffexp$hgnc_symbol),])) 
  }
  if (nrow(temp.targets) > 0) {
    measured.targets[[i]] <- temp.targets 
    # mean.effect <- plyr::ddply(temp.targets, .(Drug, Timepoint), summarize,
    #                            meanLog2FC = mean(log2FoldChange, na.rm = TRUE))
    # time.info <- rbind(time.info, mean.effect)
    # optTimes[optTimes$Drug == i,]$Timepoint <- levels(mean.effect$Timepoint)[mean.effect[which.max(abs(mean.effect$meanLog2FC)),]$Timepoint]
    # use max abs val instead of mean
    optTimes[optTimes$Drug == i,]$Timepoint <- levels(temp.targets$Timepoint)[temp.targets[which.max(abs(temp.targets$log2FoldChange)),]$Timepoint]
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
#write.csv(time.info, "meanTarget_diffexp.csv", row.names=FALSE)
write.csv(optTimes, "optimalTreatmentTimes.csv", row.names=FALSE)
write.csv(max.effect, "maximumEffect.csv", row.names=FALSE)
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
# else if no degs, then use most popular optimal time determined from measured targets
popTime <- tail(names(sort(table(optTimes$Timepoint))), 1) # source: https://stackoverflow.com/questions/12187187/how-to-retrieve-the-most-repeated-value-in-a-column-present-in-a-data-frame
opt.sens <- data.frame()
opt.sim <- data.frame()
for (i in names(drug.info)) {
  # find optimal time
  temp.time <- optTimes[optTimes$Drug == i,]$Timepoint
  drug.diffexp <- na.omit(diffexp[diffexp$Drug==i & diffexp$signif,])
  if (is.na(temp.time) & nrow(drug.diffexp) == 0) {
    temp.time <- popTime # 8h based on max log2FC
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

# make plots for each moa
for (i in unique(unlist(drug.info2))) {
  # if samples are resistant to treatment, filter for sensitivity
  # else filter for resistance if they are still sensitive to the given treatment
  
  # likewise, if drug is similar to given MOA, filter for similarity
  # else filter for opposite behavior
  
  venn.list <- list()
  strict.venn.list <- list()
  venn.colors <- list()
  temp.drugs <- names(drug.info2[drug.info2 == i])
  
  temp.sens <- opt.sens[opt.sens$DrugTreatment %in% temp.drugs,]
  given.sens <- temp.sens[temp.sens$Drug_set == i,]
  sens.dot.df <- temp.sens[temp.sens$Drug_set %in% unique(na.omit(temp.sens[temp.sens$sig,]$Drug_set)),]
  if (nrow(given.sens) > 0) {
    dir.given.sens <- mean(given.sens[given.sens$sig,]$NES, na.rm=TRUE)
    if (is.na(dir.given.sens)) {
      dir.given.sens <-  mean(given.sens$NES, na.rm=TRUE)
    } 
    # pick MOAs with opposite toxicity to that of putative MOA
    strict.sens.dot.df <- sens.dot.df[sens.dot.df$NES/dir.given.sens < 0,]
  } else if (nrow(na.omit(max.effect[max.effect$Drug %in% temp.drugs,])) > 0) {
    # else use mean of target diffexp; target upregulation suggests resistance is emerging
    dir.given.sens <- mean(max.effect[max.effect$Drug %in% temp.drugs,]$Log2FC, na.rm=TRUE)
  } else {
    dir.given.sens <- NULL
    # else pick potentially toxic MOAs if putative MOA not evaluated
    strict.sens.dot.df <- sens.dot.df[sens.dot.df$NES < 0,]
  }
  
  temp.sim <- opt.sim[opt.sim$DrugTreatment %in% temp.drugs,]
  given.sim <- temp.sim[temp.sim$Drug_set == i,]
  sim.dot.df <- temp.sim[temp.sim$Drug_set %in% unique(na.omit(temp.sim[temp.sim$sig,]$Drug_set)),]
  if (nrow(given.sim) > 0) {
    dir.given.sim <- mean(given.sim[given.sim$sig,]$NES, na.rm=TRUE)
    if (is.na(dir.given.sim)) {
      dir.given.sim <-  mean(given.sim$NES, na.rm=TRUE)
    }
    # pick MOAs with similar behavior to putative MOA
    strict.sim.dot.df <- sim.dot.df[sim.dot.df$NES/dir.given.sim > 0,]
  } else {
    dir.given.sim <- NULL
    # else pick similar MOAs if putative MOA not evaluated
    strict.sim.dot.df <- sim.dot.df[sim.dot.df$NES > 0,]
  }
  
  # pull results for each drug
  for (j in temp.drugs) {
    drug.sens <- na.omit(sens.dot.df[sens.dot.df$DrugTreatment == j,])
    drug.sim <- na.omit(sim.dot.df[sim.dot.df$DrugTreatment == j,])
    
    venn.list[[paste0("Sensitivity:\n",j)]] <- unique(na.omit(drug.sens[drug.sens$sig,]$Drug_set))
    venn.list[[paste0("Similarity:\n",j)]] <- unique(na.omit(drug.sim[drug.sim$sig,]$Drug_set))
    
    if (is.null(dir.given.sens)) {
      # default to pick potentially toxic MOAs if putative MOA not evaluated
      strict.venn.list[[paste0("Sensitivity:\n",j)]] <- unique(na.omit(drug.sens[drug.sens$sig & 
                                                                                   drug.sens$NES < 0,]$Drug_set))
      venn.colors[[paste0("Sensitivity:\n",j)]] <- "blue"
    } else {
      # else pick MOAs with opposite toxicity to that of putative MOA
      strict.venn.list[[paste0("Sensitivity:\n",j)]] <- unique(na.omit(drug.sens[drug.sens$sig & 
                                                                                   drug.sens$NES/dir.given.sens < 0,]$Drug_set))
      venn.colors[[paste0("Sensitivity:\n",j)]] <- ifelse(dir.given.sens > 0, "blue", "red")
    }
    
    if (is.null(dir.given.sim)) {
      # default to pick similar MOAs if putative MOA not evaluated
      strict.venn.list[[paste0("Similarity:\n",j)]] <- unique(na.omit(drug.sim[drug.sim$sig & 
                                                                                 drug.sim$NES > 0,]$Drug_set))
      venn.colors[[paste0("Similarity:\n",j)]] <- "red"
    } else {
      # else pick MOAs with similar behavior to putative MOA
      strict.venn.list[[paste0("Similarity:\n",j)]] <- unique(na.omit(drug.sim[drug.sim$sig & 
                                                                                 drug.sim$NES/dir.given.sim > 0,]$Drug_set)) 
      venn.colors[[paste0("Similarity:\n",j)]] <- ifelse(dir.given.sim > 0, "red", "blue")
    }
  }
  
  # venn diagrams
  if (length(temp.drugs) < 3) {
    venn.list <- venn.list[sort(names(venn.list))]
    strict.venn.list <- strict.venn.list[sort(names(strict.venn.list))]
    venn.colors <- unlist(venn.colors[names(strict.venn.list)])
    names(venn.colors) <- NULL

    ggvenn::ggvenn(venn.list, show_percentage = FALSE, set_name_size = 1, text_size = 4, show_elements=FALSE)
    ggplot2::ggsave(paste0(i,"_MOA_sensitivityANDsimilarity_vennDiagram.pdf"), width=3, height=5)
    ggvenn::ggvenn(venn.list, show_percentage = FALSE, set_name_size = 2.5, text_size = 2.5, show_elements=TRUE)
    ggplot2::ggsave(paste0(i,"_MOA_sensitivityANDsimilarity_vennDiagram_wElements.pdf"), width=5, height=5)
    
    ggvenn::ggvenn(strict.venn.list, show_percentage = FALSE, set_name_size = 1, text_size = 4, show_elements=FALSE,
                   fill_color = venn.colors)
    ggplot2::ggsave(paste0(i,"_MOA_sensitivityANDsimilarity_vennDiagram_oppSensNES_sameSimNES.pdf"), width=3, height=5)
    ggvenn::ggvenn(strict.venn.list, show_percentage = FALSE, set_name_size = 2.5, text_size = 2.5, show_elements=TRUE,
                   fill_color=venn.colors)
    ggplot2::ggsave(paste0(i,"_MOA_sensitivityANDsimilarity_vennDiagram_oppSensNES_sameSimNES_wElements.pdf"), width=5, height=5) 
  } else {
    sim.venn <- venn.list[grepl("Similarity",names(venn.list))]
    sens.venn <- venn.list[grepl("Sensitivity",names(venn.list))]
    
    strict.sim.venn <- strict.venn.list[grepl("Similarity",names(strict.venn.list))]
    strict.sens.venn <- strict.venn.list[grepl("Sensitivity",names(strict.venn.list))]
    
    sim.venn.colors <- unlist(venn.colors[names(strict.sim.venn)])
    names(sim.venn.colors) <- NULL
    sens.venn.colors <- unlist(venn.colors[names(strict.sens.venn)])
    names(sens.venn.colors) <- NULL
    
    ggvenn::ggvenn(sens.venn, show_percentage = FALSE, set_name_size = 1, text_size = 4, show_elements=FALSE)
    ggplot2::ggsave(paste0(i,"_MOA_sensitivity_vennDiagram.pdf"), width=3, height=5)
    ggvenn::ggvenn(sens.venn, show_percentage = FALSE, set_name_size = 2.5, text_size = 2.5, show_elements=TRUE)
    ggplot2::ggsave(paste0(i,"_MOA_sensitivity_vennDiagram_wElements.pdf"), width=5, height=5)
    
    ggvenn::ggvenn(sim.venn, show_percentage = FALSE, set_name_size = 1, text_size = 4, show_elements=FALSE)
    ggplot2::ggsave(paste0(i,"_MOA_similarity_vennDiagram.pdf"), width=3, height=5)
    ggvenn::ggvenn(sim.venn, show_percentage = FALSE, set_name_size = 2.5, text_size = 2.5, show_elements=TRUE)
    ggplot2::ggsave(paste0(i,"_MOA_similarity_vennDiagram_wElements.pdf"), width=5, height=5)
    
    ggvenn::ggvenn(strict.sens.venn, show_percentage = FALSE, set_name_size = 1, text_size = 4, show_elements=FALSE,
                   fill_color = sens.venn.colors)
    ggplot2::ggsave(paste0(i,"_MOA_sensitivity_vennDiagram_oppSensNES.pdf"), width=3, height=5)
    ggvenn::ggvenn(strict.sens.venn, show_percentage = FALSE, set_name_size = 2.5, text_size = 2.5, show_elements=TRUE,
                   fill_color = sens.venn.colors)
    ggplot2::ggsave(paste0(i,"_MOA_sensitivity_vennDiagram_oppSensNES_wElements.pdf"), width=5, height=5)
    
    ggvenn::ggvenn(strict.sim.venn, show_percentage = FALSE, set_name_size = 1, text_size = 4, show_elements=FALSE,
                   fill_color = sim.venn.colors)
    ggplot2::ggsave(paste0(i,"_MOA_similarity_vennDiagram_sameSimNES.pdf"), width=3, height=5)
    ggvenn::ggvenn(strict.sim.venn, show_percentage = FALSE, set_name_size = 2.5, text_size = 2.5, show_elements=TRUE,
                   fill_color = sim.venn.colors)
    ggplot2::ggsave(paste0(i,"_MOA_similarity_vennDiagram_sameSimNES_wElements.pdf"), width=5, height=5)
  }
  
  # dot plots
  simMOAs <- unique(unlist(venn.list[grepl("Similarity",names(venn.list))]))
  sensMOAs <- unique(unlist(venn.list[grepl("Sensitivity",names(venn.list))]))
  if (length(simMOAs) > 0 & length(sensMOAs) > 0) {
    if (any(sens.dot.df$minusLogFDR == "Inf")) {
      sens.dot.df[sens.dot.df$minusLogFDR == "Inf",]$minusLogFDR <- 4
    }
    sens.dot.df$Ranking <- "Sensitivity"
    sim.dot.df$minusLogFDR <- 4
    sim.dot.df[sim.dot.df$FDR_q_value !=0,]$minusLogFDR <- -log10(sim.dot.df[sim.dot.df$FDR_q_value !=0,]$FDR_q_value)
    sim.dot.df$Ranking <- "Similarity"
    keepCols <- c("MPNST", "DrugTreatment", "Drug_set", "NES", "minusLogFDR", "sig", "Ranking")
    dmea.df <- rbind(sens.dot.df[,keepCols], sim.dot.df[,keepCols])
    sigMOAs <- unique(c(simMOAs, sensMOAs))
    #sigMOAs <- unique(c(sens.dot.df[sens.dot.df$sig,]$Drug_set, sim.dot.df[sim.dot.df$sig,]$Drug_set)) # 46
    sharedMOAs <- sens.dot.df$Drug_set[sens.dot.df$Drug_set %in% sim.dot.df$Drug_set] # 52
    sigSharedMOAs <- sigMOAs[sigMOAs %in% sharedMOAs]
    if (length(sigSharedMOAs) > 0) {
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
          x = DrugTreatment, y = Drug_set, color = NES,
          size = minusLogFDR
        )
      ) + facet_grid(MPNST ~ Ranking) +
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
      ggplot2::ggsave(paste0(i,"_Compiled_DMEA_dotPlot_sensitivityANDsimilarity.pdf"), width=4.2, height=4)
      
      if (all(grepl(" inhibitor", sigOrder))) {
        dot.df$Inhibitor <- sub(" inhibitor", "", dot.df$Drug_set)
        sigOrder2 <- sub(" inhibitor", "", sigOrder)
        longNames <- c("Aurora kinase", "Bromodomain", "Sterol demethylase")
        names(longNames) <- c("AURK", "BET", "CYP51")
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
        ) + facet_grid(MPNST ~ Ranking) +
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
        ggplot2::ggsave(paste0(i,"_Compiled_DMEA_dotPlot_sensitivityANDsimilarity_v2.pdf"), width=3.4, height=4)
      } else {
        ggplot2::ggsave(paste0(i,"_Compiled_DMEA_dotPlot_sensitivityANDsimilarity_wider.pdf"), width=5, height=4)
        ggplot2::ggsave(paste0(i,"_Compiled_DMEA_dotPlot_sensitivityANDsimilarity_widerTaller.pdf"), width=5, height=5)
      }
    }
    
    # dot plots of potentially synergistic MOAs
    if (any(strict.sens.dot.df$minusLogFDR == "Inf")) {
      strict.sens.dot.df[strict.sens.dot.df$minusLogFDR == "Inf",]$minusLogFDR <- 4 
    }
    strict.sens.dot.df$Ranking <- "Sensitivity"
    strict.sim.dot.df$minusLogFDR <- 4
    strict.sim.dot.df[strict.sim.dot.df$FDR_q_value !=0,]$minusLogFDR <- 
      -log10(strict.sim.dot.df[strict.sim.dot.df$FDR_q_value !=0,]$FDR_q_value)
    strict.sim.dot.df$Ranking <- "Similarity"
    keepCols <- c("MPNST", "DrugTreatment", "Drug_set", "NES", "minusLogFDR", "sig", "Ranking")
    strict.dmea.df <- rbind(strict.sens.dot.df[,keepCols], strict.sim.dot.df[,keepCols])
    simMOAs <- unique(unlist(strict.venn.list[grepl("Similarity",names(strict.venn.list))]))
    sensMOAs <- unique(unlist(strict.venn.list[grepl("Sensitivity",names(strict.venn.list))]))
    synMOAs <- simMOAs[simMOAs %in% sensMOAs]
    if (length(synMOAs) > 0) {
      mean.dmea.df <- plyr::ddply(strict.dmea.df, .(Drug_set), summarize,
                                  absNES = mean(abs(NES), na.rm=TRUE))
      sigOrder <- mean.dmea.df[order(mean.dmea.df$absNES),]$Drug_set
      sigOrder <- sigOrder[sigOrder %in% synMOAs]
      dot.df <- dmea.df[dmea.df$Drug_set %in% synMOAs,]
      library(ggplot2)
      maxAbsNES <- max(abs(dot.df$NES))
      dot.plot <- ggplot2::ggplot(
        dot.df,
        ggplot2::aes(
          x = DrugTreatment, y = Drug_set, color = NES,
          size = minusLogFDR
        )
      ) + facet_grid(MPNST ~ Ranking) +
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
      ggplot2::ggsave(paste0(i,"_Compiled_DMEA_dotPlot_sensitivityANDsimilarity_potentialSynergy.pdf"), width=4.2, height=4)
      
      if (all(grepl(" inhibitor", sigOrder))) {
        dot.df$Inhibitor <- sub(" inhibitor", "", dot.df$Drug_set)
        sigOrder2 <- sub(" inhibitor", "", sigOrder)
        longNames <- c("Aurora kinase", "Bromodomain", "Sterol demethylase")
        names(longNames) <- c("AURK", "BET", "CYP51")
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
        ) + facet_grid(MPNST ~ Ranking) +
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
        ggplot2::ggsave(paste0(i,"_Compiled_DMEA_dotPlot_sensitivityANDsimilarity_potentialSynergy_v2.pdf"), width=3.4, height=4)
      } else {
        ggplot2::ggsave(paste0(i,"_Compiled_DMEA_dotPlot_sensitivityANDsimilarity_potentialSynergy_wider.pdf"), width=5, height=4)
        ggplot2::ggsave(paste0(i,"_Compiled_DMEA_dotPlot_sensitivityANDsimilarity_potentialSynergy_widerTaller.pdf"), width=5, height=5)
      }
    } 
  }
}
