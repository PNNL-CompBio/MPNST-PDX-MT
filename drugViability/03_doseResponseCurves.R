# analyzing IncuCyte results
library(plyr); library(dplyr); library(tidyr);library(ggplot2);library(synapser)
library(drc) # curve source: answer by greenjune: https://stackoverflow.com/questions/36780357/plotting-dose-response-curves-with-ggplot2-and-drc
# Sara shared these 2 links: https://stackoverflow.com/questions/68209998/plot-drc-using-ggplot; https://forum.posit.co/t/extract-out-points-in-dose-response-curve-produced-by-drc/159433

synapser::synLogin()
setwd("~/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/MPNST-PDX-MT/drugViability")
dataPath <- "~/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/MPNST-PDX-MT/drugViability"

#### prep data ####
rel.conf <- read.csv("mpnst_combo_drug_response.csv") # relative viability
musyc.scores <- read.csv("Deconvolved_musyc_20250812.csv") # musyc synergy
bliss <- read.csv("bliss_results.csv") # bliss synergy

# calculate mean and sd for each drug & dose combo, mpnst, time combo; also preserve sample and chr8q columns
mean.conf <- plyr::ddply(rel.conf, .(sample, drug1, drug1.conc, drug2, drug2.conc), summarize,
                         meanGROWTH = mean(effect, na.rm=TRUE),
                         sdGROWTH = sd(effect, na.rm=TRUE))
mean.conf$time <- sub(".*_","",mean.conf$sample)
mean.conf$time <- sub("hours","",mean.conf$time)
mean.conf$time <- as.numeric(mean.conf$time)/24

# fix sample names
full.names <- c("MN-4","JH-2-079c", "JH-2-002", "WU-225","WU-356", "MN-2")
short.names <- c("N-4","79c", "002", "225", "356", "N-2")
for (i in 1:length(short.names)) {
  mean.conf[grepl(short.names[[i]],mean.conf$sample),]$sample <- full.names[[i]]
}
prc2.loss <- c("JH-2-079c", "JH-2-002", "WU-225", "WU-365", "MN-2")
prc2.wt <- c("MN-4")
mean.conf$PRC2 <- "Deficient"
mean.conf[mean.conf$sample %in% prc2.wt,]$PRC2 <- "Intact"

mean.conf$conc1 <- mean.conf$drug1.conc
mean.conf$conc2 <- mean.conf$drug2.conc
bliss.tested <- merge(bliss, mean.conf, by=c("sample","time","drug1","drug2","conc1","conc2")) # 2529 down from 10287 bliss
max.bliss.tested <- plyr::ddply(bliss.tested, .(drugCombo, sample, time), summarize,
                         maxBliss = max(Bliss_synergy)) # 111
write.csv(max.bliss.tested,"maxBlissTested.csv", row.names=FALSE)

bliss.musyc.tested <- merge(max.bliss.tested, musyc.scores, by=c("drugCombo","sample","time"))
results <- merge(mean.conf, musyc.scores, by=c("sample","time","drug1","drug2"), all.x=TRUE)
results <- merge(results, max.bliss.tested, by=c("sample","time","drugCombo"), all.x=TRUE)
write.csv(results,"results_viabilityBlissMusyc.csv", row.names=FALSE)
results$timeD <- paste0(results$time,"d")

##### dose response curves #####
results$conc12.ratio <- signif(as.numeric(results$drug1.conc/results$drug2.conc),3)
results$drugCombo <- paste0(results$drug1,"+",results$drug2)
rel.conf$drugCombo <- paste0(rel.conf$drug1,"+",rel.conf$drug2)

combos <- unique(results$drugCombo) # 21
for (c in combos) {
  combo.res <- results[results$drugCombo == c,]
  d1 <- unique(combo.res$drug1)
  d2 <- unique(combo.res$drug2)
  if (nrow(combo.res) > 0) {
    # generate plots
    temp.cols2 <- c("blue","black","red")
    if (length(unique(combo.res$conc12.ratio)) == 4) {temp.cols2 <- c("blue","black","grey","red")}
    mid.color2 <- "black"
    
    conf.plot <- ggplot(combo.res, aes(x=conc, y=100*meanGROWTH, color=as.factor(conc12.ratio))) +
      facet_grid(timeD ~ sample) +
      geom_smooth(linetype="dashed", se=FALSE, method=drc::drm, method.args=list(fct=L.4(), se=FALSE)) +
      scale_x_continuous(transform="log10") + geom_point() + 
      geom_errorbar(aes(ymin=100*(meanGROWTH-sdGROWTH), ymax=100*(meanGROWTH+sdGROWTH))) +
      ggtitle(paste0(d1," + ",d2))+scale_color_manual(values=temp.cols2) +
      theme_classic(base_size=12) + labs(x="Concentration (uM)", y = "% Viability",
                                         color = paste("Ratio of",d1,"\nto",d2)) +
      theme(plot.title=element_text(hjust=0.5,face="bold"), 
            axis.text.x=element_text(angle=45, hjust=1, vjust=1)) + 
      geom_text(data=combo.res, mapping=aes(x=Inf, y=Inf, fontface="plain",
                                            label=paste0("Log|",as.character("\u03b1|=("),
                                                         signif(log_alpha12,2),",",signif(log_alpha21,2),
                                                         "),\nMax Bliss=",
                                                         signif(maxBliss,2))),
                                        hjust=1, vjust=1, show.legend=FALSE, 
                                        colour=mid.color2)
    ggsave(paste0(d2,"_",d1,"_viability_log10_ratios_MuSyCmaxBlissTested_p_blueBlackRed_noR2plain.svg"),conf.plot,width=12,height=9, device="svg") # was height 4
    ggsave(paste0(d2,"_",d1,"_viability_log10_ratios_MuSyCmaxBlissTested_p_blueBlackRed_noR2plain.png"),conf.plot,width=12,height=9, device="png") # was height 4
    ggsave(paste0(d2,"_",d1,"_viability_log10_ratios_MuSyCmaxBlissTested_p_blueBlackRed_noR2plain_h4w12.svg"),conf.plot,width=12,height=4, device="svg") # was height 4
    ggsave(paste0(d2,"_",d1,"_viability_log10_ratios_MuSyCmaxBlissTested_p_blueBlackRed_noR2plain_h4w12.png"),conf.plot,width=12,height=4, device="png") # was height 4
    ggsave(paste0(d2,"_",d1,"_viability_log10_ratios_MuSyCmaxBlissTested_p_blueBlackRed_noR2plain_h4w10.svg"),conf.plot,width=10,height=4, device="svg") # was height 4
    ggsave(paste0(d2,"_",d1,"_viability_log10_ratios_MuSyCmaxBlissTested_p_blueBlackRed_noR2plain_h4w10.png"),conf.plot,width=10,height=4, device="png") # was height 4
    ggsave(paste0(d2,"_",d1,"_viability_log10_ratios_MuSyCmaxBlissTested_p_blueBlackRed_noR2plain_h4w5.svg"),conf.plot,width=5,height=4, device="svg") # was height 4
    ggsave(paste0(d2,"_",d1,"_viability_log10_ratios_MuSyCmaxBlissTested_p_blueBlackRed_noR2plain_h4w5.png"),conf.plot,width=5,height=4, device="png") # was height 4
  }
}

#### single agent curves ####
sing.conf <- read.table("mpnst_drug_response.tsv", sep="\t", header=TRUE)
sing.conf$X <- FALSE
library(drc) # curve source: answer by greenjune: https://stackoverflow.com/questions/36780357/plotting-dose-response-curves-with-ggplot2-and-drc
# Sara shared these 2 links: https://stackoverflow.com/questions/68209998/plot-drc-using-ggplot; https://forum.posit.co/t/extract-out-points-in-dose-response-curve-produced-by-drc/159433
mean.sing <- plyr::ddply(sing.conf, .(DOSE, improve_sample_id, Drug, time), summarize,
                         meanGROWTH=mean(GROWTH),
                         sdGROWTH=sd(GROWTH))
mean.sing$timeD <- paste0(mean.sing$time/24,"d")
conf.plot <- ggplot(mean.sing, aes(x=DOSE, y=meanGROWTH, color=improve_sample_id)) +
  facet_grid(timeD ~ Drug) +
  geom_smooth(linetype="dashed", se=FALSE, method=drc::drm, method.args=list(fct=L.4(), se=FALSE)) +
  scale_x_continuous(transform="log10") + geom_point() + 
  geom_errorbar(aes(ymin=(meanGROWTH-sdGROWTH), ymax=(meanGROWTH+sdGROWTH))) +
  #ggtitle(paste0(d1," + ",d2))+scale_color_manual(values=temp.cols2) +
  theme_classic(base_size=12) + labs(x="Concentration (uM)", y = "% Viability", color="MPNST") +
  theme(plot.title=element_text(hjust=0.5,face="bold"), axis.text.x=element_text(angle=45, hjust=1, vjust=1))
conf.plot # something is off- looks like MN2 isn't relative viability, same for some MN4
ggsave("relativeViability_singleAgents.pdf", conf.plot, width=18, height=3)

top.sing <- mean.sing[mean.sing$Drug %in% c("vorinostat","mirdametinib","RMC4630"),]
top.sing$Drug <- factor(top.sing$Drug, levels=c("vorinostat","mirdametinib","RMC4630"))
conf.plot <- ggplot(top.sing, 
                    aes(x=DOSE, y=meanGROWTH, color=improve_sample_id)) +
  facet_grid(Drug ~ timeD) +
  geom_smooth(linetype="dashed", se=FALSE, method=drc::drm, method.args=list(fct=L.4(), se=FALSE)) +
  scale_x_continuous(transform="log10") + geom_point() + 
  geom_errorbar(aes(ymin=(meanGROWTH-sdGROWTH), ymax=(meanGROWTH+sdGROWTH))) +
  #ggtitle(paste0(d1," + ",d2))+scale_color_manual(values=temp.cols2) +
  theme_classic(base_size=12) + labs(x="Concentration (uM)", y = "% Viability", color="MPNST") +
  theme(plot.title=element_text(hjust=0.5,face="bold"), axis.text.x=element_text(angle=45, hjust=1, vjust=1))
conf.plot # something is off- looks like MN2 isn't relative viability, same for some MN4
ggsave("relativeViability_top3Drugs.pdf", conf.plot, width=6, height=6)
ggsave("relativeViability_top3Drugs_w4h4.pdf", conf.plot, width=4, height=4)
ggsave("relativeViability_top3Drugs_w5h5.pdf", conf.plot, width=5, height=5)