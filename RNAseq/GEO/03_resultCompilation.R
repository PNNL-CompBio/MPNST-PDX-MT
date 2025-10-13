# MEK+HDACi result compilation
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

#### compare DMEA results with similarity or sensitivity ####
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
sens <- read.csv("Compiled/DMEA_results.csv")
sim <- read.csv("DMEA_Cmap_L1000/DMEA_Cmap_L1000.csv")
sim$sig <- FALSE
sim[sim$p_value <= 0.05 & sim$FDR_q_value <= 0.25,]$sig <- TRUE
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
