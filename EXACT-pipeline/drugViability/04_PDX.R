# compile DMEA for potential synergy
rm(list=ls())
library(plyr);library(dplyr);library(synapser);library(data.table);library(ggplot2);library(patchwork)
synapser::synLogin()

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
dot.df <- na.omit(tumor.size)
dot.df$SD <- as.numeric(dot.df$SD)
ggplot(dot.df,aes(x=`Treatment.Day..`,y=`Mean.Size`, color=Treatment)) +
  geom_point() + geom_errorbar(aes(ymin=`Mean.Size`-SD, ymax=`Mean.Size`+SD))+
  geom_smooth(se=FALSE, linetype="dashed")+
  theme_classic() + labs(y="Mean Tumor Size (mm^3)", x="Treatment Duration (Days)") +
  scale_color_manual(values=c("black","red","blue","forestgreen"), breaks=c("Vehicle","Mirdametinib","Vorinostat","Mirdametinib + Vorinostat")) +
  facet_wrap(.~MPNST, ncol=2)
#ggsave("mirdaVorinostat_PDX.pdf", width=6, height=3) # was width 4
write.csv(tumor.size,"PDX_mirdametinibVorinostat_meanTumorSize.csv", row.names=FALSE)
tumor.size <- read.csv("PDX_mirdametinibVorinostat_meanTumorSize.csv")

dot.df <- na.omit(tumor.size[tumor.size$MPNST %in% c("MN-2","WU-225"),])
dot.df$SD <- as.numeric(dot.df$SD)
ggplot(dot.df,
       aes(x=`Treatment.Day..`,y=`Mean.Size`, color=Treatment)) +
  geom_point() + geom_errorbar(aes(ymin=`Mean.Size`-SD, ymax=`Mean.Size`+SD))+
  geom_smooth(se=FALSE, linetype="dashed")+
  theme_classic() + labs(y=paste0("Mean Tumor Size (m",expression("m^3"),")"), x="Treatment Duration (Days)") +
  scale_color_manual(values=c("grey","red","blue","black"),
                     breaks=c("Vehicle","Mirdametinib","Vorinostat","Mirdametinib + Vorinostat")) +
  facet_wrap(.~MPNST, ncol=2)
#ggsave("mirdaVorinostat_PDX_MN2-WU225.pdf", width=6, height=3) # was width 4
ggsave("fig_7b_mirdaVorinostat_PDX_MN2-WU225_wide.pdf", width=12, height=3) # was width 4
#
# combo.t <- t.test(tumor.size[tumor.size$Treatment %in% c("Mirdametinib","Vorinostat"),]$`Mean.Size`,
#                   tumor.size[tumor.size$Treatment == "Mirdametinib + Vorinostat",]$`Mean.Size`, "greater")
# combo.t # p = 0.07263
#
# mek.combo.t <- t.test(tumor.size[tumor.size$Treatment == "Mirdametinib",]$`Mean.Size`,
#                   tumor.size[tumor.size$Treatment == "Mirdametinib + Vorinostat",]$`Mean.Size`, "greater", paired=TRUE)
# mek.combo.t # p = 0.1716
#
# vor.combo.t <- t.test(tumor.size[tumor.size$Treatment == "Vorinostat",]$`Mean.Size`,
#                       tumor.size[tumor.size$Treatment == "Mirdametinib + Vorinostat",]$`Mean.Size`, "greater", paired=TRUE)
# vor.combo.t # error due to unequal vector lengths
#
# veh.mek.t <- t.test(tumor.size[tumor.size$Treatment == "Vehicle",]$`Mean.Size`,
#                       tumor.size[tumor.size$Treatment == "Mirdametinib",]$`Mean.Size`, "greater", paired=TRUE)
# veh.mek.t # error due to unequal vector lengths
#
# veh.vor.t <- t.test(tumor.size[tumor.size$Treatment == "Vehicle",]$`Mean.Size`,
#                       tumor.size[tumor.size$Treatment == "Vorinostat",]$`Mean.Size`, "greater", paired=TRUE)
# veh.vor.t # error due to unequal vector lengths

mpnsts <- unique(tumor.size$MPNST)
p.df <- data.frame()
for (m in mpnsts) {
  t.df <- data.frame(test=c("singleGreaterThanCombo","mirdaGreaterThanCombo",
                            "vorinGreaterThanCombo","vehicleGreaterThanMirda",
                            "vehicleGreaterThanVorin","vehicleGreaterThanCombo"),p=NA)
  m.tumor.size <- tumor.size[tumor.size$MPNST==m,]
  combo.t <- t.test(m.tumor.size[m.tumor.size$Treatment %in% c("Mirdametinib","Vorinostat"),]$`Mean.Size`,
                    m.tumor.size[m.tumor.size$Treatment == "Mirdametinib + Vorinostat",]$`Mean.Size`, "greater")
  combo.t # p = 0.07263
  t.df[t.df$test=="singleGreaterThanCombo",]$p <- combo.t$p.value

  mek.combo.t <- t.test(m.tumor.size[m.tumor.size$Treatment == "Mirdametinib",]$`Mean.Size`,
                        m.tumor.size[m.tumor.size$Treatment == "Mirdametinib + Vorinostat",]$`Mean.Size`, "greater")
  mek.combo.t # p = 0.1716
  t.df[t.df$test=="mirdaGreaterThanCombo",]$p <- mek.combo.t$p.value

  vor.combo.t <- t.test(m.tumor.size[m.tumor.size$Treatment == "Vorinostat",]$`Mean.Size`,
                        m.tumor.size[m.tumor.size$Treatment == "Mirdametinib + Vorinostat",]$`Mean.Size`, "greater")
  vor.combo.t # error due to unequal vector lengths
  t.df[t.df$test=="vorinGreaterThanCombo",]$p <- vor.combo.t$p.value

  veh.mek.t <- t.test(m.tumor.size[m.tumor.size$Treatment == "Vehicle",]$`Mean.Size`,
                      m.tumor.size[m.tumor.size$Treatment == "Mirdametinib",]$`Mean.Size`, "greater")
  veh.mek.t # error due to unequal vector lengths
  t.df[t.df$test=="vehicleGreaterThanMirda",]$p <- veh.mek.t$p.value

  veh.vor.t <- t.test(m.tumor.size[m.tumor.size$Treatment == "Vehicle",]$`Mean.Size`,
                      m.tumor.size[m.tumor.size$Treatment == "Vorinostat",]$`Mean.Size`, "greater")
  veh.vor.t # error due to unequal vector lengths
  t.df[t.df$test=="vehicleGreaterThanVorin",]$p <- veh.vor.t$p.value

  veh.combo.t <- t.test(m.tumor.size[m.tumor.size$Treatment == "Vehicle",]$`Mean.Size`,
                      m.tumor.size[m.tumor.size$Treatment == "Mirdametinib + Vorinostat",]$`Mean.Size`, "greater")
  veh.combo.t # error due to unequal vector lengths
  t.df[t.df$test=="vehicleGreaterThanCombo",]$p <- veh.combo.t$p.value
  t.df$MPNST <- m
  p.df <- rbind(p.df, t.df)
}
p.df$q <- qvalue::qvalue(p.df$p, pi0=1)$qvalues
write.csv(p.df,"PDX_mirdaVorin_pValues.csv",row.names=FALSE)

p.df.h <- p.df[p.df$MPNST!="JH-2-002",]
p.df.h$q <- qvalue::qvalue(p.df.h$p, pi0=1)$qvalues
write.csv(p.df.h,"PDX_mirdaVorin_pValues_HirbeLab.csv",row.names=FALSE)
#
