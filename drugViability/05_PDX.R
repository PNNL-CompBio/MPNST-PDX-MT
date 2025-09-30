# compile DMEA for potential synergy
rm(list=ls())
library(plyr);library(dplyr);library(synapser);library(data.table);library(ggplot2);library(patchwork)
synapser::synLogin()
setwd("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/MPNST-PDX-MT/drugViability")
dir.create("PDX")
setwd("PDX")

#### PDX data: MEKi+HDACi ####
#tumor.size <- read.csv("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/MPNST-PDX-MT/Wu225 Mirda + Vorinostat.csv")
wu225 <- read.csv(synapser::synGet("syn68900596")$path)
wu225[wu225$Treatment=="Mirda",]$Treatment <- "Mirdametinib"
wu225[wu225$Treatment=="Mirda + Vorinostat",]$Treatment <- "Mirdametinib + Vorinostat"
mn2 <- read.csv(synapser::synGet("syn69953578")$path)
colnames(mn2)[2] <- "Treatment"
wu356 <- read.csv(synapser::synGet("syn69953579")$path)
colnames(wu356)[2] <- "Treatment"
jh2002 <- read.csv(synapser::synGet("syn69953662")$path)
colnames(jh2002)[2] <- "Treatment"
colnames(wu225) <- colnames(mn2)
tumor.size <- na.omit(rbind(wu225, mn2, wu356, jh2002))

# get colors for mirda and vorinostat
library(RColorBrewer); library(scales)
cols1=c("#000000",brewer.pal(8,'Dark2'),brewer.pal(15-8,'Set2'),"mediumorchid1", "cornflowerblue", "#004B4B", "#4B0026")
drug.info <- list("palbociclib" = "CDK inhibitor", "ribociclib" = "CDK inhibitor",
                  "trabectedin" = "Chemotherapy", "Ifosfamide" = "DNA alkylating agent",
                  "decitabine" = "DNMT inhibitor", "vorinostat" = "HDAC inhibitor",
                  "mirdametinib" = "MEK inhibitor", "selumetinib" = "MEK inhibitor",
                  "trametinib" = "MEK inhibitor", "capmatinib" = "MET inhibitor",
                  "olaparib" = "PARP inhibitor", "RMC4630" = "SHP2 inhibitor",
                  "TNO155" = "SHP2 inhibitor", "doxorubicin" = "TOP inhibitor",
                  "irinotecan" = "TOP inhibitor", "verteporfin" = "YAP inhibitor",
                  "pexidartinib" = "KIT inhibitor", "IAG933" = "YAP inhibitor", "SN38" = "TOP inhibitor")
names(cols1) <- c("DMSO", names(drug.info))
scales::show_col(cols1[c("DMSO","mirdametinib","vorinostat")]) # black, brown, yellow?
ggplot(na.omit(tumor.size),aes(x=`Treatment.Day..`,y=`Mean.Size`, color=Treatment)) + 
  geom_point() + geom_errorbar(aes(ymin=`Mean.Size`-SD, ymax=`Mean.Size`+SD))+
  geom_smooth(se=FALSE, linetype="dashed")+
  theme_classic() + labs(y="Mean Tumor Size (mm^3)", x="Treatment Duration (Days)") +
  scale_color_manual(values=c("black","red","blue","forestgreen"), breaks=c("Vehicle","Mirdametinib","Vorinostat","Mirdametinib + Vorinostat")) + 
  facet_wrap(.~MPNST, ncol=2)
setwd("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/MPNST-PDX-MT")
ggsave("mirdaVorinostat_PDX.pdf", width=6, height=3) # was width 4
write.csv(tumor.size,"PDX_mirdametinibVorinostat_meanTumorSize.csv", row.names=FALSE)

combo.t <- t.test(tumor.size[tumor.size$Treatment %in% c("Mirda","Vorinostat"),]$`Mean.Size`,
                  tumor.size[tumor.size$Treatment == "Mirda + Vorinostat",]$`Mean.Size`, "greater")
combo.t # p = 0.000649

mek.combo.t <- t.test(tumor.size[tumor.size$Treatment == "Mirda",]$`Mean.Size`,
                  tumor.size[tumor.size$Treatment == "Mirda + Vorinostat",]$`Mean.Size`, "greater", paired=TRUE)
mek.combo.t # p = 0.001507

vor.combo.t <- t.test(tumor.size[tumor.size$Treatment == "Vorinostat",]$`Mean.Size`,
                      tumor.size[tumor.size$Treatment == "Mirda + Vorinostat",]$`Mean.Size`, "greater", paired=TRUE)
vor.combo.t # p = 0.002504

veh.mek.t <- t.test(tumor.size[tumor.size$Treatment == "Vehicle",]$`Mean.Size`,
                      tumor.size[tumor.size$Treatment == "Mirda",]$`Mean.Size`, "greater", paired=TRUE)
veh.mek.t # p = 0.006131

veh.vor.t <- t.test(tumor.size[tumor.size$Treatment == "Vehicle",]$`Mean.Size`,
                      tumor.size[tumor.size$Treatment == "Vorinostat",]$`Mean.Size`, "greater", paired=TRUE)
veh.vor.t # p = 0.004812

#### mirda combos predicted ####
mirda.sens <- sens[sens$DrugTreatment=="Mirdametinib" & sens$sig & sens$NES < 0,]
mirda.mean <- plyr::ddply(mirda.sens, .(Drug_set), summarize,
                          meanNES=mean(NES),
                          medianNES=median(NES),
                          sdNES=sd(NES),
                          N=dplyr::n(),
                          N_MPNST = length(unique(MPNST)),
                          MPNST = paste0(sort(unique(MPNST)), collapse=", "))
mirda.mean$Tested <- "plain"
mirda.mean[mirda.mean$Drug_set %in% drug.info2,]$Tested <- "bold"
ggplot(mirda.mean, aes(y=reorder(Drug_set, -meanNES), x=meanNES, fill=MPNST)) + 
  geom_bar(stat="identity") + theme_classic() + scale_y_discrete(position="right") +
  theme(axis.title.y=element_blank(), legend.position="bottom") + labs(x="Mean NES") +
  scale_fill_manual(values=RColorBrewer::brewer.pal(length(unique(mirda.mean$MPNST)),"Dark2"), breaks=unique(mirda.mean$MPNST))
  #scale_fill_manual(values=c("black","darkgrey","lightgrey"), breaks=c("JH-2-002 & MN-2","JH-2-002","MN-2"))
ggsave("Mirda_DMEA_predictions.pdf", width=5, height=5)

ggplot(mirda.mean, aes(x=reorder(Drug_set, meanNES), y=meanNES, fill=MPNST)) + 
  geom_bar(stat="identity") + theme_classic() + 
  theme(axis.title.x=element_blank(), legend.position="top", axis.text.x=element_text(angle=45, vjust=1, hjust=1)) + labs(y="Mean NES") +
  scale_fill_manual(values=RColorBrewer::brewer.pal(length(unique(mirda.mean$MPNST)),"Dark2"), breaks=unique(mirda.mean$MPNST))
ggsave("Mirda_DMEA_predictions_horizontal.pdf", width=5, height=5)

mirda.mean[mirda.mean$Drug_set == "protein tyrosine kinase inhibitor",]$Drug_set <- "Protein TKI"
mirda.mean[mirda.mean$Drug_set == "tyrosine kinase inhibitor",]$Drug_set <- "TKI"
mirda.mean[grepl("PDGFR", mirda.mean$Drug_set),]$Drug_set <- "PDGFRi"
mirda.mean[grepl("Aurora", mirda.mean$Drug_set),]$Drug_set <- "AURKi"
mirda.mean[grepl("bromodomain", mirda.mean$Drug_set),]$Drug_set <- "BRDi"
mirda.mean[grepl("adenosine", mirda.mean$Drug_set),]$Drug_set <- "P1i"
mirda.mean[grepl("glutamate", mirda.mean$Drug_set),]$Drug_set <- "GluRi"
mirda.mean[grepl("src", mirda.mean$Drug_set),]$Drug_set <- "SRCi"
mirda.mean[grepl("retinoid ", mirda.mean$Drug_set),]$Drug_set <- "RR Agonist"
mirda.mean[grepl("inhibitor", mirda.mean$Drug_set),]$Drug_set <- sub(" inhibitor","i",mirda.mean[grepl("inhibitor", mirda.mean$Drug_set),]$Drug_set)
ggplot(mirda.mean, aes(x=reorder(Drug_set, meanNES), y=meanNES, fill=MPNST)) + 
  geom_bar(stat="identity") + theme_classic() + 
  theme(axis.title.x=element_blank(), legend.position="top", axis.text.x=element_text(angle=45, vjust=1, hjust=1)) + labs(y="Mean NES") +
  scale_fill_manual(values=RColorBrewer::brewer.pal(length(unique(mirda.mean$MPNST)),"Dark2"), breaks=unique(mirda.mean$MPNST))
ggsave("Mirda_DMEA_predictions_horizontal_shortNames.pdf", width=5, height=3)
ggsave("Mirda_DMEA_predictions_horizontal_shortNames_lessWide.pdf", width=4, height=3)

ggplot(mirda.mean, aes(y=reorder(Drug_set, -meanNES), x=meanNES, fill=MPNST)) + 
  geom_bar(stat="identity") + theme_classic() + scale_y_discrete(position="right") +
  theme(axis.title.y=element_blank(), legend.position="bottom") + labs(x="Mean NES") +
  scale_fill_manual(values=RColorBrewer::brewer.pal(length(unique(mirda.mean$MPNST)),"Dark2"), breaks=unique(mirda.mean$MPNST))
ggsave("Mirda_DMEA_predictions_shortNames.pdf", width=4, height=5)
ggsave("Mirda_DMEA_predictions_shortNames_smaller.pdf", width=3, height=4)
