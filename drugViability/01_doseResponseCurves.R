# analyzing IncuCyte results
library(plyr); library(dplyr); library(tidyr);library(ggplot2);library(colorspace);library(synapser)
synapser::synLogin()
setwd("~/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/MPNST-PDX-MT/drugViability")
dataPath <- "~/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/MPNST-PDX-MT/drugViability"

rel.conf <- read.csv("mpnst_combo_drug_response.csv")
musyc.scores <- read.csv("musyc_20250812.csv")

# try looking at Bliss scores
bliss <- data.frame()

#dir.create("bliss")
blissFiles <- synapser::as.list(synapser::synGetChildren('syn68639935', list("file"), sortBy = 'NAME'))
for (i in 1:length(blissFiles)){
  test.bliss <- read.csv(synapser::synGet(blissFiles[[i]]$id)$path)
  bliss <- rbind(bliss, test.bliss)
}
bliss$drugCombo <- paste0(bliss$drug1, "+", bliss$drug2)
bliss$sample <- bliss$PDX
bliss$time <- bliss$timePoint..hr./24
write.csv(bliss, "bliss_results.csv", row.names=FALSE)
bliss <- read.csv("bliss_results.csv")

# calculate mean and sd for each drug & dose combo, mpnst, time combo; also preserve sample and chr8q columns
mean.conf <- plyr::ddply(rel.conf, .(sample, drug1, drug1.conc, drug2, drug2.conc), summarize,
                         meanGROWTH = mean(effect, na.rm=TRUE),
                         sdGROWTH = sd(effect, na.rm=TRUE))
mean.conf$time <- sub(".*_","",mean.conf$sample)
mean.conf$time <- sub("hours","",mean.conf$time)
mean.conf$time <- as.numeric(mean.conf$time)/24

musyc.scores$time <- sub(".*_","",musyc.scores$sample)
musyc.scores$time <- sub("hours","",musyc.scores$time)
musyc.scores$time <- as.numeric(musyc.scores$time)/24

# fix sample names
full.names <- c("MN-4","JH-2-079c", "JH-2-002", "WU-225","WU-356", "MN-2")
short.names <- c("N-4","79c", "002", "225", "356", "N-2")
for (i in 1:length(short.names)) {
  mean.conf[grepl(short.names[[i]],mean.conf$sample),]$sample <- full.names[[i]]
  musyc.scores[grepl(short.names[[i]],musyc.scores$sample),]$sample <- full.names[[i]]
}
prc2.loss <- c("JH-2-079c", "JH-2-002", "WU-225", "WU-365", "MN-2")
prc2.wt <- c("MN-4")
mean.conf$PRC2 <- "Deficient"
mean.conf[mean.conf$sample %in% prc2.wt,]$PRC2 <- "Intact"

# restore full drug names for musyc data
drugs <- unique(c(mean.conf$drug1, mean.conf$drug2))
short.drugs <- substr(drugs, 2, 4)
colnames(musyc.scores)[2:3] <- c("drug1","drug2")
for (i in 1:length(drugs)) {
  if (any(musyc.scores$drug1 == short.drugs[[i]])) {
    musyc.scores[musyc.scores$drug1 == short.drugs[[i]],]$drug1 <- drugs[[i]] 
  }
  if (any(musyc.scores$drug2 == short.drugs[[i]])) {
    musyc.scores[musyc.scores$drug2 == short.drugs[[i]],]$drug2 <- drugs[[i]] 
  }
}
write.csv(musyc.scores, "Deconvolved_musyc_20250812.csv", row.names=FALSE) # duplicated and changed name to musyc.csv for synapse
musyc.scores <- read.csv("Deconvolved_musyc_20250812.csv")

# make scatter plot of a12 vs a21
dir.create(paste0("curves_",Sys.Date()))
setwd(paste0("curves_",Sys.Date()))
musyc.scores$drugCombo <- paste0(musyc.scores$drug1,"+",musyc.scores$drug2)
ggplot(musyc.scores, aes(x=log_alpha21, y=log_alpha12, shape=sample, color=drugCombo, size=R2)) +
  geom_point()+theme_minimal() + 
  #ggrepel::geom_text_repel(aes(label=paste0(time,"d"))) + 
  labs(x="Drug 2's Effect on Potency of Drug 1", y="Drug 2's Effect on Potency of Drug 1")

ggplot(musyc.scores[musyc.scores$log_alpha12>0 & musyc.scores$log_alpha21>0,], 
       aes(x=log_alpha21, y=log_alpha12, shape=sample, color=drugCombo, size=R2)) +
  geom_point()+theme_minimal() + 
  #ggrepel::geom_text_repel(aes(label=paste0(time,"d"))) + 
  labs(x="Drug 2's Effect on Potency of Drug 1", y="Drug 2's Effect on Potency of Drug 1")

musyc.scores$shortCombo <- paste0(tolower(substr(musyc.scores$drug1, 1,4)), "+",
                                  tolower(substr(musyc.scores$drug2, 1,4)))
ggplot(musyc.scores[musyc.scores$log_alpha12>0 & musyc.scores$log_alpha21>0,], 
       aes(x=log_alpha21, y=log_alpha12, shape=as.factor(time), color=drugCombo, size=R2)) +
  geom_point()+theme_classic() + facet_wrap(.~sample)+ scale_size_continuous(breaks=c(0.3,0.5,0.7)) +
  ggrepel::geom_text_repel(aes(label=shortCombo), box.padding = 0.5) + 
  labs(x="Drug 2's Effect on Potency of Drug 1", y="Drug 2's Effect on Potency of Drug 1", shape="Time (d)")
ggsave("musyc_potentCombos_v2.pdf",width=12, height=9)
#ggsave("musyc_potentCombos_classicTheme.pdf",width=12, height=6)

musyc.scores$meanLogAlpha <- NA
for (i in 1:nrow(musyc.scores)) {
  musyc.scores$meanLogAlpha[i] <- mean(c(musyc.scores$log_alpha12[i], musyc.scores$log_alpha21[i])) 
}
maxLogAlpha <- max(musyc.scores$meanLogAlpha)
musyc.scores$drugCombo <- paste0(musyc.scores$drug1,'+',musyc.scores$drug2)
ggplot(musyc.scores, aes(x=sample, y=reorder(drugCombo, meanLogAlpha), fill=meanLogAlpha)) + geom_tile(stat="identity") + 
  facet_wrap(.~paste0(time,"d"))+theme_classic() + scale_fill_gradient(low="grey",high="red",limits=c(0,8)) + 
  theme(axis.title = element_blank(), axis.text.x = element_text(angle=45, vjust=1, hjust=1)) + labs(fill="Mean\nLog|Alpha|")
ggsave("musyc_meanLogAlpha_heatmap_positive.pdf", width=4,height=4)

ggplot(musyc.scores, aes(x=sample, y=reorder(drugCombo, meanLogAlpha), fill=meanLogAlpha)) + geom_tile(stat="identity") + 
  facet_wrap(.~paste0(time,"d"))+theme_classic() + scale_fill_gradient2(low="blue",mid="grey",high="red",limits=c(-8,8)) + 
  theme(axis.title = element_blank(), axis.text.x = element_text(angle=45, vjust=1, hjust=1)) + labs(fill="Mean\nLog|Alpha|")
ggsave("musyc_meanLogAlpha_heatmap.pdf", width=5,height=4)


bliss.musyc <- merge(bliss, musyc.scores, by=c("drugCombo","sample","time"))
max.bliss <- plyr::ddply(bliss, .(drugCombo, sample, time), summarize,
                         maxBliss = max(Bliss_synergy)) # 125
write.csv(max.bliss,"maxBliss.csv", row.names=FALSE)

mean.conf$conc1 <- mean.conf$drug1.conc
mean.conf$conc2 <- mean.conf$drug2.conc
bliss.tested <- merge(bliss, mean.conf, by=c("sample","time","drug1","drug2","conc1","conc2")) # 2529 down from 10287 bliss
max.bliss.tested <- plyr::ddply(bliss.tested, .(drugCombo, sample, time), summarize,
                         maxBliss = max(Bliss_synergy)) # 111
write.csv(max.bliss.tested,"maxBlissTested.csv", row.names=FALSE)

bliss.musyc <- merge(max.bliss, musyc.scores, by=c("drugCombo","sample","time"))
ggplot(bliss.musyc, aes(x=meanLogAlpha, y=maxBliss, shape=sample, color=drugCombo)) + #, size=R2)) +
  geom_point()+theme_minimal() + 
  #ggrepel::geom_text_repel(aes(label=paste0(time,"d"))) + 
  labs(x="Mean MuSyC Potency", y="Max Bliss Synergy")
ggsave("musyc_meanLogAlpha_vs_max_bliss_synergy.pdf", width=10, height=8)

ggplot(bliss.musyc, aes(x=meanLogAlpha, y=maxBliss, shape=sample, color=drugCombo)) + #, size=R2)) +
  geom_point()+theme_minimal() + facet_grid(paste0(time,"d")~sample)+
  #ggrepel::geom_text_repel(aes(label=paste0(time,"d"))) + 
  labs(x="Mean MuSyC Potency", y="Max Bliss Synergy")
ggsave("musyc_meanLogAlpha_vs_max_bliss_synergy_Faceted.pdf", width=12, height=7)

cor.test(bliss.musyc$meanLogAlpha, bliss.musyc$maxBliss) # r = 0.018, p=0.8247
cor.test(bliss.musyc$meanLogAlpha, bliss.musyc$maxBliss, method="spearman") # rho = 0.043, p = 0.5939

bliss.musyc.tested <- merge(max.bliss.tested, musyc.scores, by=c("drugCombo","sample","time"))
ggplot(bliss.musyc.tested, aes(x=meanLogAlpha, y=maxBliss, shape=sample, color=drugCombo)) + #, size=R2)) +
  geom_point()+theme_minimal() + 
  #ggrepel::geom_text_repel(aes(label=paste0(time,"d"))) + 
  labs(x="Mean MuSyC Potency", y="Max Bliss Synergy")
ggsave("musyc_meanLogAlpha_vs_max_bliss_synergy_tested.pdf", width=10, height=8)

ggplot(bliss.musyc.tested, aes(x=meanLogAlpha, y=maxBliss, shape=sample, color=drugCombo)) + #, size=R2)) +
  geom_point()+theme_minimal() + facet_grid(paste0(time,"d")~sample)+
  #ggrepel::geom_text_repel(aes(label=paste0(time,"d"))) + 
  labs(x="Mean MuSyC Potency", y="Max Bliss Synergy")
ggsave("musyc_meanLogAlpha_vs_max_bliss_synergy_tested_Faceted.pdf", width=12, height=7)

cor.test(bliss.musyc.tested$meanLogAlpha, bliss.musyc.tested$maxBliss) # r=0.06, p=0.218
cor.test(bliss.musyc.tested$meanLogAlpha, bliss.musyc.tested$maxBliss, method="spearman") # rho = 0.118, p=0.1475

results <- merge(mean.conf, musyc.scores, by=c("sample","time","drug1","drug2"), all.x=TRUE)
results <- merge(results, max.bliss.tested, by=c("sample","time","drugCombo"), all.x=TRUE)
write.csv(results,"results_viabilityBlissMusyc.csv", row.names=FALSE)
#results <- read.csv(file.path(dataPath,"curves_2025-08-04/results_viabilityBlissMusyc.csv"))
library(drc) # curve source: answer by greenjune: https://stackoverflow.com/questions/36780357/plotting-dose-response-curves-with-ggplot2-and-drc
# Sara shared these 2 links: https://stackoverflow.com/questions/68209998/plot-drc-using-ggplot; https://forum.posit.co/t/extract-out-points-in-dose-response-curve-produced-by-drc/159433
results$timeD <- paste0(results$time,"d")
library(RColorBrewer)
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
scales::show_col(cols1[drugs])
scales::show_col(cols1)
drugs[!(drugs %in% names(cols1))] # none missing
cols <- colorspace::darken(cols1, 0.3)
scales::show_col(cols)

##### ratio centric #####
dir.create(paste0("curves_",Sys.Date()))
setwd(paste0("curves_",Sys.Date()))
#results$combo.conc <- NULL
results$conc12.ratio <- signif(as.numeric(results$drug1.conc/results$drug2.conc),3)
results$drugCombo <- paste0(results$drug1,"+",results$drug2)
library(plyr);library(dplyr);library(tidyr)
results$potentialSynergyMuSyC <- FALSE
results[!is.na(results$meanLogAlpha) & results$log_alpha12>0 & 
          results$log_alpha21>0,]$potentialSynergyMuSyC <- TRUE
results$potentialSynergyBliss <- FALSE
results[!is.na(results$maxBliss) & results$maxBliss>0,]$potentialSynergyBliss <- TRUE
results$textFace <- "plain"
results[results$potentialSynergyMuSyC | results$potentialSynergyBliss,]$textFace <- "bold"

rel.conf$drugCombo <- paste0(rel.conf$drug1,"+",rel.conf$drug2)

all.p.df <- data.frame()
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
    
    # t-test: if viability lower with combo vs. single agents?
    combo.conf <- rel.conf[rel.conf$drugCombo == c,]
    sample <- c(unique(combo.conf$sample),"Overall")
    p.df <- data.frame(drugCombo = c, sample, comboMean = NA, singleMean = NA, 
                       N_combo = NA, N_single = NA, p = NA, t = NA)
    combo.df <- combo.conf[combo.conf$drug1.conc!=0 &
                             combo.conf$drug2.conc!=0,]
    sing.df <- combo.conf[combo.conf$drug1.conc==0 |
                             combo.conf$drug2.conc==0 & 
                            !(combo.conf$drug1.conc==0 &
                                combo.conf$drug2.conc==0),]
    
    if (nrow(combo.df) > 0 & nrow(sing.df) > 0) {
      combo.test <- t.test(combo.df$effect, sing.df$effect, "less")
      p.df[p.df$sample=="Overall",]$comboMean <- combo.test$estimate[["mean of x"]]
      p.df[p.df$sample=="Overall",]$singleMean <- combo.test$estimate[["mean of y"]]
      p.df[p.df$sample=="Overall",]$N_combo <- length(combo.df$effect)
      p.df[p.df$sample=="Overall",]$N_single <- length(sing.df$effect)
      p.df[p.df$sample=="Overall",]$p <- combo.test$p.value
      p.df[p.df$sample=="Overall",]$t <- combo.test$statistic[["t"]]
      for (s in unique(combo.conf$sample)) {
        scombo.df <- combo.df[combo.df$sample == s,]
        ssing.df <- sing.df[sing.df$sample == s,]
        if (nrow(scombo.df) > 0 & nrow(ssing.df) > 0) {
          combo.test <- t.test(scombo.df$effect, ssing.df$effect, "less") 
          p.df[p.df$sample==s,]$comboMean <- combo.test$estimate[["mean of x"]]
          p.df[p.df$sample==s,]$singleMean <- combo.test$estimate[["mean of y"]]
          p.df[p.df$sample==s,]$N_combo <- length(scombo.df$effect)
          p.df[p.df$sample==s,]$N_single <- length(ssing.df$effect)
          p.df[p.df$sample==s,]$p <- combo.test$p.value
          p.df[p.df$sample==s,]$t <- combo.test$statistic[["t"]]
        }
      }
      all.p.df <- rbind(all.p.df, p.df)
      p.df <- p.df[p.df$sample != "Overall",] %>% 
        tidyr::separate_wider_delim(sample, "_", names=c("sample", "time"))
      p.df$time <- as.numeric(sub("hours","",p.df$time))/24
      
      combo.resp <- merge(combo.res, p.df, by=c("drugCombo","sample","time"), all.x=TRUE)
      combo.resp <- combo.resp[,c("drugCombo","sample","timeD","drug1","drug2",
                                  "drug1.conc","drug2.conc","meanGROWTH",
                                  "sdGROWTH","log_alpha12","log_alpha21","R2",
                                  "maxBliss","p","textFace","conc12.ratio")]
      drug1.df <- combo.resp[combo.resp$drug2.conc==0,]
      drug1.df$drug2.conc <- NULL
      drug1.df$drug2 <- NULL
      drug1.df$drugCombo <- NULL
      colnames(drug1.df) <- c("sample","timeD","drug", "conc","meanGROWTH",
                              "sdGROWTH","log_alpha12","log_alpha21","R2",
                              "maxBliss","p","textFace","conc12.ratio")
      
      drug2.df <- combo.resp[combo.resp$drug1.conc==0,]
      drug2.df$drug1.conc <- NULL
      drug2.df$drug1 <- NULL
      drug2.df$drugCombo <- NULL
      colnames(drug2.df) <- c("sample","timeD","drug", "conc","meanGROWTH",
                              "sdGROWTH","log_alpha12","log_alpha21","R2",
                              "maxBliss","p","textFace","conc12.ratio")
      
      combo.resp2 <- rbind(drug1.df, drug2.df)
      combo.ratios <- unique(combo.resp$conc12.ratio[combo.resp$conc12.ratio != 0 &
                                                        combo.resp$conc12.ratio != "Inf"])
      for (combo.ratio in combo.ratios) {
        drug12.df <- combo.resp[combo.resp$conc12.ratio == combo.ratio,]
        if (mean(drug1.df$meanGROWTH) < mean(drug2.df$meanGROWTH)) {
          drug12.df$combo.conc <- drug12.df$drug1.conc
        } else {
          drug12.df$combo.conc <- drug12.df$drug2.conc
        }
        drug12.df$drug1.conc <- NULL
        drug12.df$drug1 <- NULL
        drug12.df$drug2.conc <- NULL
        drug12.df$drug2 <- NULL
        colnames(drug12.df) <- c("drug","sample","timeD","meanGROWTH",
                                 "sdGROWTH","log_alpha12","log_alpha21","R2",
                                 "maxBliss","p","textFace","conc12.ratio","conc")
        combo.resp2 <- rbind(combo.resp2, drug12.df[,colnames(combo.resp2)])
      }
    } else {
      combo.resp2 <- data.frame()
    }
    
    if (nrow(combo.resp2) > 0) {

      conf.plot <- ggplot(combo.resp2, aes(x=conc, y=100*meanGROWTH, color=as.factor(conc12.ratio))) +
        facet_grid(timeD ~ sample) +
        geom_smooth(linetype="dashed", se=FALSE, method=drc::drm, method.args=list(fct=L.4(), se=FALSE)) +
        scale_x_continuous(transform="log10") + geom_point() + 
        geom_errorbar(aes(ymin=100*(meanGROWTH-sdGROWTH), ymax=100*(meanGROWTH+sdGROWTH))) +
        ggtitle(paste0(d1," + ",d2))+scale_color_manual(values=temp.cols2) +
        theme_classic(base_size=12) + labs(x="Concentration (uM)", y = "% Viability",
                                           color = paste("Ratio of",d1,"\nto",d2)) +
        theme(plot.title=element_text(hjust=0.5,face="bold"), 
              axis.text.x=element_text(angle=45, hjust=1, vjust=1))

      conf.plotC <- conf.plot + geom_text(data=combo.resp2, 
                                          mapping=aes(x=Inf, y=Inf, fontface="plain",
                                                      label=paste0("Log|",as.character("\u03b1|=("),
                                                                   signif(log_alpha12,2),",",signif(log_alpha21,2),
                                                                   "),\nMax Bliss=",
                                                                   signif(maxBliss,2))),
                                          hjust=1, vjust=1, show.legend=FALSE, 
                                          colour=mid.color2)
      ggsave(paste0(d2,"_",d1,"_viability_log10_ratios_MuSyCmaxBlissTested_p_blueBlackRed_noR2plain.svg"),conf.plotC,width=12,height=9, device="svg") # was height 4
      ggsave(paste0(d2,"_",d1,"_viability_log10_ratios_MuSyCmaxBlissTested_p_blueBlackRed_noR2plain.png"),conf.plotC,width=12,height=9, device="png") # was height 4
      ggsave(paste0(d2,"_",d1,"_viability_log10_ratios_MuSyCmaxBlissTested_p_blueBlackRed_noR2plain_h4w12.svg"),conf.plotC,width=12,height=4, device="svg") # was height 4
      ggsave(paste0(d2,"_",d1,"_viability_log10_ratios_MuSyCmaxBlissTested_p_blueBlackRed_noR2plain_h4w12.png"),conf.plotC,width=12,height=4, device="png") # was height 4
      ggsave(paste0(d2,"_",d1,"_viability_log10_ratios_MuSyCmaxBlissTested_p_blueBlackRed_noR2plain_h4w10.svg"),conf.plotC,width=10,height=4, device="svg") # was height 4
      ggsave(paste0(d2,"_",d1,"_viability_log10_ratios_MuSyCmaxBlissTested_p_blueBlackRed_noR2plain_h4w10.png"),conf.plotC,width=10,height=4, device="png") # was height 4
      ggsave(paste0(d2,"_",d1,"_viability_log10_ratios_MuSyCmaxBlissTested_p_blueBlackRed_noR2plain_h4w5.svg"),conf.plotC,width=5,height=4, device="svg") # was height 4
      ggsave(paste0(d2,"_",d1,"_viability_log10_ratios_MuSyCmaxBlissTested_p_blueBlackRed_noR2plain_h4w5.png"),conf.plotC,width=5,height=4, device="png") # was height 4
    }
  }
}
write.csv(all.p.df, "tTests_viability_combo_lessThan_single.csv", row.names=FALSE)

#### compare musyc alpha to DMEA ####
##### drug corr #####
synapser::synLogin()
drug.corr <- read.csv(synapser::synGet("syn65672266")$path)
drug.corr[drug.corr$DrugTreatment == "Trabectidin",]$DrugTreatment <- "Trabectedin"

shared.drugs <- drugs[tolower(drugs) %in% tolower(c(drug.corr$Drug, drug.corr$DrugTreatment))] # 9: mirdametinib, olaparib, TNO155, capmatinib, doxorubicin, palbociclib, vorinostat, irinotecan, trabectedin
shared.combos <- tidyr::expand_grid(shared.drugs, shared.drugs)
shared.combos$drugCombo <- paste0(shared.combos$shared.drugs...1,"+",shared.combos$shared.drugs...2)
shared.combos$drugCombo2 <- paste0(shared.combos$shared.drugs...2,"+",shared.combos$shared.drugs...1)
shared.combos <- unique(c(shared.combos$drugCombo, shared.combos$drugCombo2))
shared.combos <- shared.combos[shared.combos %in% combos] # 23 / 27 tested combos

shared.corr <- drug.corr[tolower(drug.corr$Drug) %in% tolower(drugs) & 
                           tolower(drug.corr$DrugTreatment) %in% tolower(drugs),]
shared.corr$drugCombo <- NA
for (i in 1:nrow(shared.corr)) {
  d1 <- shared.corr$DrugTreatment[i]
  d2 <- shared.corr$Drug[i]
  
  # make lowercase unless TNO155, IAG933, or SN38
  d1 <- ifelse(d1 %in% c("TNO155","IAG933","SN38"),d1,tolower(d1))
  d2 <- ifelse(d2 %in% c("TNO155","IAG933","SN38"),d2,tolower(d2))
  
  # find right combo order (d1+d2 or d2+d1)
  temp.combo <- c(paste0(d1,"+",d2),paste0(d2,"+",d1))
  temp.combo <- temp.combo[temp.combo %in% shared.combos]
  if (length(temp.combo == 1)) {shared.corr$drugCombo[i] <- temp.combo} else if (length(temp.combo) > 1){
    warning("more than 1 combo for",d1, "and", d2)
  }
}
colnames(shared.corr)[2] <- "sample"
shared.corr <- na.omit(shared.corr[,c("sample","drugCombo",
                                      "Pearson.est","Pearson.p","Pearson.q",
                                      "Spearman.est","Spearman.p","Spearman.q","N")])
drug.corr.res <- merge(shared.corr, results, by=c("drugCombo","sample")) # no time overlap

# note: only mirdametinib+doxorubicin is significantly correlated in MN-2 (Pearson.q<0.05) which has Pearson.est=0.1915
drug.corr.res <- na.omit(drug.corr.res)
write.csv(drug.corr.res, "drug_dmea_vs_synergy.csv", row.names = FALSE)
# create colors using color in-between drug colors
drug.corr.res$color <- NA
for (i in shared.combos) {
  if (i %in% drug.corr.res$drugCombo) {
    ds <- strsplit(i,"[+]")[[1]]
    col1 <- cols[[ds[1]]]
    col2 <- cols[[ds[2]]]
    temp.cols <- colorRampPalette(c(col1, col2))(3)
    drug.corr.res[drug.corr.res$drugCombo == i,]$color <- temp.cols[2] 
  }
}
combo.cols <- drug.corr.res$color
names(combo.cols) <- drug.corr.res$drugCombo
combo.cols <- unique(combo.cols) # 9
drug.corr.res <- dplyr::distinct(drug.corr.res)

syn.scores <- c("Maximum Bliss Synergy Score" = "maxBliss", "Mean MuSyC Potency Score" = "meanLogAlpha")
for (i in names(syn.scores)) {
  drug.corr.res$rank <- drug.corr.res[,syn.scores[i]]
  # signed
  musyc.drug.corr <- cor.test(drug.corr.res$rank, drug.corr.res$Pearson.est)
  musyc.drug.corr$estimate # 0.2382369 
  musyc.drug.corr$p.value # 1.293807e-12
  musyc.drug.corrsp <- cor.test(drug.corr.res$rank, drug.corr.res$Pearson.est,method="spearman")
  musyc.drug.corrsp$estimate # 0.2414861
  musyc.drug.corrsp$p.value # 6.260826e-13
  
  musyc.drug.corr2 <- cor.test(drug.corr.res$rank, drug.corr.res$Spearman.est)
  musyc.drug.corr2$estimate # 0.3405962 
  musyc.drug.corr2$p.value # 6.573243e-25
  musyc.drug.corr2sp <- cor.test(drug.corr.res$rank, drug.corr.res$Spearman.est,method="spearman")
  musyc.drug.corr2sp$estimate # 0.3188854 
  musyc.drug.corr2sp$p.value # 7.132333e-22
  
  if (musyc.drug.corr$p.value <= 0.05) {
    stats_pearson <- substitute(
      r == est * "," ~ ~"p" ~ "=" ~ p,
      list(
        est = format(as.numeric(musyc.drug.corr$estimate), digits = 2),
        p = format(musyc.drug.corr$p.value, digits = 2)
      )
    )
    ggplot(drug.corr.res, aes(x=rank, y=Pearson.est)) + # could set shape=sample but only sample is JH-2-002 
      geom_point(aes(color=drugCombo), show.legend =TRUE) + theme_minimal() + scale_color_manual(values=combo.cols, breaks=names(combo.cols),labels=names(combo.cols)) + 
      geom_smooth(method="lm",linetype="dashed",color="black") + ggplot2::geom_text(
        x = -Inf, y = Inf, vjust = "inward", hjust = "inward",
        colour = "black", parse = TRUE, 
        label = as.character(as.expression(stats_pearson)), size = 8
      ) +
      labs(x=i,
           y="DMEA Sensitivity Prediction")
    ggsave(paste0("drugResults_",i,"_vs_dmea_PearsonEst.pdf"),width=5,height=5)
  }
  
  if (musyc.drug.corrsp$p.value <= 0.05) {
    stats_spearman <- substitute(
      rho == est * "," ~ ~"p" ~ "=" ~ p,
      list(
        est = format(as.numeric(musyc.drug.corrsp$estimate), digits = 2),
        p = format(musyc.drug.corrsp$p.value, digits = 2)
      )
    )
    ggplot(drug.corr.res, aes(x=rank, y=Pearson.est)) + # could set shape=sample but only sample is JH-2-002 
      geom_point(aes(color=drugCombo), show.legend =TRUE) + theme_minimal() + scale_color_manual(values=combo.cols, breaks=names(combo.cols),labels=names(combo.cols)) + 
      geom_smooth(method="lm",linetype="dashed",color="black") + ggplot2::geom_text(
        x = -Inf, y = Inf, vjust = "inward", hjust = "inward",
        colour = "black", parse = TRUE, 
        label = as.character(as.expression(stats_spearman)), size = 8
      ) +
      labs(x=i,
           y="DMEA Sensitivity Prediction")
    ggsave(paste0("drugResults_",i,"_vs_dmea_PearsonEst_spearman.pdf"),width=5,height=5)
  }
  
  if (musyc.drug.corr2$p.value <= 0.05) {
    stats_pearson <- substitute(
      r == est * "," ~ ~"p" ~ "=" ~ p,
      list(
        est = format(as.numeric(musyc.drug.corr2$estimate), digits = 2),
        p = format(musyc.drug.corr2$p.value, digits = 2)
      )
    )
    ggplot(drug.corr.res, aes(x=rank, y=Spearman.est)) + # could set shape=sample but only sample is JH-2-002 
      geom_point(aes(color=drugCombo), show.legend =TRUE) + theme_minimal() + scale_color_manual(values=combo.cols, breaks=names(combo.cols),labels=names(combo.cols)) + 
      geom_smooth(method="lm",linetype="dashed",color="black") + ggplot2::geom_text(
        x = -Inf, y = Inf, vjust = "inward", hjust = "inward",
        colour = "black", parse = TRUE, 
        label = as.character(as.expression(stats_pearson)), size = 8
      ) +
      labs(x=i,
           y="DMEA Sensitivity Prediction")
    ggsave(paste0("drugResults_",i,"_vs_dmea_SpearmanEst.pdf"),width=5,height=5)
  }
  
  if (musyc.drug.corr2sp$p.value <= 0.05) {
    stats_spearman <- substitute(
      rho == est * "," ~ ~"p" ~ "=" ~ p,
      list(
        est = format(as.numeric(musyc.drug.corr2sp$estimate), digits = 2),
        p = format(musyc.drug.corr2sp$p.value, digits = 2)
      )
    )
    ggplot(drug.corr.res, aes(x=rank, y=Spearman.est)) + # could set shape=sample but only sample is JH-2-002 
      geom_point(aes(color=drugCombo), show.legend =TRUE) + theme_minimal() + scale_color_manual(values=combo.cols, breaks=names(combo.cols),labels=names(combo.cols)) + 
      geom_smooth(method="lm",linetype="dashed",color="black") + ggplot2::geom_text(
        x = -Inf, y = Inf, vjust = "inward", hjust = "inward",
        colour = "black", parse = TRUE, 
        label = as.character(as.expression(stats_spearman)), size = 8
      ) +
      labs(x=i,
           y="DMEA Sensitivity Prediction")
    ggsave(paste0("drugResults_",i,"_vs_dmea_SpearmanEst_spearman.pdf"),width=5,height=5)
  }
  
  # abs drug corr
  musyc.drug.corr <- cor.test(drug.corr.res$rank, abs(drug.corr.res$Pearson.est))
  musyc.drug.corr$estimate # 0.2382369 
  musyc.drug.corr$p.value # 1.293807e-12
  musyc.drug.corrsp <- cor.test(drug.corr.res$rank, abs(drug.corr.res$Pearson.est),method="spearman")
  musyc.drug.corrsp$estimate # 0.2414861
  musyc.drug.corrsp$p.value # 6.260826e-13
  
  musyc.drug.corr2 <- cor.test(drug.corr.res$rank, abs(drug.corr.res$Spearman.est))
  musyc.drug.corr2$estimate # 0.3405962 
  musyc.drug.corr2$p.value # 6.573243e-25
  musyc.drug.corr2sp <- cor.test(drug.corr.res$rank, abs(drug.corr.res$Spearman.est),method="spearman")
  musyc.drug.corr2sp$estimate # 0.3188854 
  musyc.drug.corr2sp$p.value # 7.132333e-22
  
  if (musyc.drug.corr$p.value <= 0.05) {
    stats_pearson <- substitute(
      r == est * "," ~ ~"p" ~ "=" ~ p,
      list(
        est = format(as.numeric(musyc.drug.corr$estimate), digits = 2),
        p = format(musyc.drug.corr$p.value, digits = 2)
      )
    )
    ggplot(drug.corr.res, aes(x=rank, y=abs(Pearson.est))) + # could set shape=sample but only sample is JH-2-002 
      geom_point(aes(color=drugCombo), show.legend =TRUE) + theme_minimal() + scale_color_manual(values=combo.cols, breaks=names(combo.cols),labels=names(combo.cols)) + 
      geom_smooth(method="lm",linetype="dashed",color="black") + ggplot2::geom_text(
        x = -Inf, y = Inf, vjust = "inward", hjust = "inward",
        colour = "black", parse = TRUE, 
        label = as.character(as.expression(stats_pearson)), size = 8
      ) +
      labs(x=i,
           y="Absolute DMEA Sensitivity Prediction")
    ggsave(paste0("drugResults_",i,"_vs_dmea_absPearsonEst.pdf"),width=5,height=5)
  }
  
  if (musyc.drug.corrsp$p.value <= 0.05) {
    stats_spearman <- substitute(
      rho == est * "," ~ ~"p" ~ "=" ~ p,
      list(
        est = format(as.numeric(musyc.drug.corrsp$estimate), digits = 2),
        p = format(musyc.drug.corrsp$p.value, digits = 2)
      )
    )
    ggplot(drug.corr.res, aes(x=rank, y=abs(Pearson.est))) + # could set shape=sample but only sample is JH-2-002 
      geom_point(aes(color=drugCombo), show.legend =TRUE) + theme_minimal() + scale_color_manual(values=combo.cols, breaks=names(combo.cols),labels=names(combo.cols)) + 
      geom_smooth(method="lm",linetype="dashed",color="black") + ggplot2::geom_text(
        x = -Inf, y = Inf, vjust = "inward", hjust = "inward",
        colour = "black", parse = TRUE, 
        label = as.character(as.expression(stats_spearman)), size = 8
      ) +
      labs(x=i,
           y="Absolute DMEA Sensitivity Prediction")
    ggsave(paste0("drugResults_",i,"_vs_dmea_absPearsonEst_spearman.pdf"),width=5,height=5)
  }
  
  if (musyc.drug.corr2$p.value <= 0.05) {
    stats_pearson <- substitute(
      r == est * "," ~ ~"p" ~ "=" ~ p,
      list(
        est = format(as.numeric(musyc.drug.corr2$estimate), digits = 2),
        p = format(musyc.drug.corr2$p.value, digits = 2)
      )
    )
    ggplot(drug.corr.res, aes(x=rank, y=abs(Spearman.est))) + # could set shape=sample but only sample is JH-2-002 
      geom_point(aes(color=drugCombo), show.legend =TRUE) + theme_minimal() + scale_color_manual(values=combo.cols, breaks=names(combo.cols),labels=names(combo.cols)) + 
      geom_smooth(method="lm",linetype="dashed",color="black") + ggplot2::geom_text(
        x = -Inf, y = Inf, vjust = "inward", hjust = "inward",
        colour = "black", parse = TRUE, 
        label = as.character(as.expression(stats_pearson)), size = 8
      ) +
      labs(x=i,
           y="Absolute DMEA Sensitivity Prediction")
    ggsave(paste0("drugResults_",i,"_vs_dmea_absSpearmanEst.pdf"),width=5,height=5)
  }
  
  if (musyc.drug.corr2sp$p.value <= 0.05) {
    stats_spearman <- substitute(
      rho == est * "," ~ ~"p" ~ "=" ~ p,
      list(
        est = format(as.numeric(musyc.drug.corr2sp$estimate), digits = 2),
        p = format(musyc.drug.corr2sp$p.value, digits = 2)
      )
    )
    ggplot(drug.corr.res, aes(x=rank, y=abs(Spearman.est))) + # could set shape=sample but only sample is JH-2-002 
      geom_point(aes(color=drugCombo), show.legend =TRUE) + theme_minimal() + scale_color_manual(values=combo.cols, breaks=names(combo.cols),labels=names(combo.cols)) + 
      geom_smooth(method="lm",linetype="dashed",color="black") + ggplot2::geom_text(
        x = -Inf, y = Inf, vjust = "inward", hjust = "inward",
        colour = "black", parse = TRUE, 
        label = as.character(as.expression(stats_spearman)), size = 8
      ) +
      labs(x=i,
           y="Absolute DMEA Sensitivity Prediction")
    ggsave(paste0("drugResults_",i,"_vs_dmea_absSpearmanEst_spearman.pdf"),width=5,height=5) 
  }
}

##### moas #####
moa.df <- read.csv(synapser::synGet("syn65672258")$path)
moa.df[moa.df$DrugTreatment=="Trabectidin",]$DrugTreatment <- "Trabectedin"
shared.moas <- drug.info[tolower(drug.info) %in% tolower(moa.df$Drug_set)] # 10
# are CSF1R inhibitors evaluated? (pexidartinib also inhibits CSF1R in addition to KIT)
any(grepl("CSF1R", moa.df$Drug_set, ignore.case=TRUE)) # FALSE

# moa.combos <- c()
# for (c in shared.combos2) {
#   ds <- strsplit(c,"[+]")[[1]]
#   moa1 <- drug.info[[ds[1]]]
#   moa2 <- drug.info[[ds[2]]]
#   moa.combos <- c(moa.combos, paste0(moa1,"+",moa2))
# }

shared.moa.df <- moa.df[tolower(moa.df$DrugTreatment) %in% tolower(drugs) & 
                          tolower(moa.df$Drug_set) %in% tolower(shared.moas) & 
                          moa.df$MPNST %in% results$sample,]
shared.moa.df$moa <- NA
for (d in shared.drugs) {
  shared.moa.df[tolower(shared.moa.df$DrugTreatment) == tolower(d),]$moa <- drug.info[[d]]
}
length(unique(shared.moa.df$moa)) # 8
shared.moas2 <- shared.moas[shared.moas %in% shared.moa.df$moa] # 8

results2 <- results[,c("sample","drug1","drug2","meanLogAlpha","maxBliss","R2","timeD")]
results2[,c("moa1","moa2","moaCombo")] <- NA
for (d in drugs) {
  if (any(results2$drug1 == d)) {
    results2[results$drug1==d,]$moa1 <- drug.info[[d]]
  }
  if (any(results2$drug2 == d)) {
    results2[results$drug2==d,]$moa2 <- drug.info[[d]]
  }
}

results3 <- results2[results2$moa1 %in% shared.moas2 &
                       results2$moa2 %in% shared.moas2,]
results3$moaCombo <- paste0(results3$moa1,"+",results3$moa2)
#results3 <- results3[results3$moaCombo %in% moa.combos,]
moa.combos <- unique(results3$moaCombo) # 7

shared.moa.df$moaCombo <- NA
for (c in moa.combos) {
  cmoas <- strsplit(c, "[+]")[[1]]
  shared.moa.df[shared.moa.df$moa %in% cmoas & shared.moa.df$Drug_set %in% cmoas &
                  shared.moa.df$moa != shared.moa.df$Drug_set,]$moaCombo <- c
}
shared.moa.df2 <- na.omit(shared.moa.df)
colnames(shared.moa.df2)[2] <- "sample"
shared.moa.res <- merge(shared.moa.df2, results3, by=c("sample","moaCombo"))
shared.moa.res$moaComboShort <- gsub(" inhibitor","i",shared.moa.res$moaCombo)
shared.moa.res <- dplyr::distinct(shared.moa.res) # 40 rows
write.csv(shared.moa.res, "moa_dmea_vs_synergy.csv", row.names = FALSE)

syn.scores <- c("Maximum Bliss Synergy Score" = "maxBliss", "Mean MuSyC Potency Score" = "meanLogAlpha")
for (i in names(syn.scores)) {
  shared.moa.res$rank <- shared.moa.res[,syn.scores[i]]
  # signed
  musyc.drug.corr <- cor.test(shared.moa.res$rank, shared.moa.res$NES)
  musyc.drug.corr$estimate # 0.2382369 
  musyc.drug.corr$p.value # 1.293807e-12
  musyc.drug.corrsp <- cor.test(shared.moa.res$rank, shared.moa.res$NES,method="spearman")
  musyc.drug.corrsp$estimate # 0.2414861
  musyc.drug.corrsp$p.value # 6.260826e-13
  
  stats_pearson <- substitute(
    r == est * "," ~ ~"p" ~ "=" ~ p,
    list(
      est = format(as.numeric(musyc.drug.corr$estimate), digits = 2),
      p = format(musyc.drug.corr$p.value, digits = 2)
    )
  )
  ggplot(shared.moa.res, aes(x=rank, y=NES)) + # could set shape=sample but only sample is JH-2-002 
    geom_point(aes(color=moaCombo), show.legend =TRUE) + theme_minimal() + scale_color_manual(values=combo.cols, breaks=names(combo.cols),labels=names(combo.cols)) + 
    geom_smooth(method="lm",linetype="dashed",color="black") + ggplot2::geom_text(
      x = -Inf, y = Inf, vjust = "inward", hjust = "inward",
      colour = "black", parse = TRUE, 
      label = as.character(as.expression(stats_pearson)), size = 8
    ) +
    labs(x=i,
         y="DMEA Sensitivity Prediction")
  ggsave(paste0("moaResults_",i,"_vs_dmea_NES.pdf"),width=5,height=5)
  
  stats_spearman <- substitute(
    rho == est * "," ~ ~"p" ~ "=" ~ p,
    list(
      est = format(as.numeric(musyc.drug.corrsp$estimate), digits = 2),
      p = format(musyc.drug.corrsp$p.value, digits = 2)
    )
  )
  ggplot(shared.moa.res, aes(x=rank, y=NES)) + # could set shape=sample but only sample is JH-2-002 
    geom_point(aes(color=moaCombo), show.legend =TRUE) + theme_minimal() + scale_color_manual(values=combo.cols, breaks=names(combo.cols),labels=names(combo.cols)) + 
    geom_smooth(method="lm",linetype="dashed",color="black") + ggplot2::geom_text(
      x = -Inf, y = Inf, vjust = "inward", hjust = "inward",
      colour = "black", parse = TRUE, 
      label = as.character(as.expression(stats_spearman)), size = 8
    ) +
    labs(x=i,
         y="DMEA Sensitivity Prediction")
  ggsave(paste0("moaResults_",i,"_vs_dmea_NES_spearman.pdf"),width=5,height=5)
  
  # abs drug corr
  musyc.drug.corr <- cor.test(shared.moa.res$rank, abs(shared.moa.res$NES))
  musyc.drug.corr$estimate # 0.2382369 
  musyc.drug.corr$p.value # 1.293807e-12
  musyc.drug.corrsp <- cor.test(shared.moa.res$rank, abs(shared.moa.res$NES),method="spearman")
  musyc.drug.corrsp$estimate # 0.2414861
  musyc.drug.corrsp$p.value # 6.260826e-13
  
  stats_pearson <- substitute(
    r == est * "," ~ ~"p" ~ "=" ~ p,
    list(
      est = format(as.numeric(musyc.drug.corr$estimate), digits = 2),
      p = format(musyc.drug.corr$p.value, digits = 2)
    )
  )
  ggplot(shared.moa.res, aes(x=rank, y=abs(NES))) + # could set shape=sample but only sample is JH-2-002 
    geom_point(aes(color=moaCombo), show.legend =TRUE) + theme_minimal() + scale_color_manual(values=combo.cols, breaks=names(combo.cols),labels=names(combo.cols)) + 
    geom_smooth(method="lm",linetype="dashed",color="black") + ggplot2::geom_text(
      x = -Inf, y = Inf, vjust = "inward", hjust = "inward",
      colour = "black", parse = TRUE, 
      label = as.character(as.expression(stats_pearson)), size = 8
    ) +
    labs(x=i,
         y="Absolute DMEA Sensitivity Prediction")
  ggsave(paste0("moaResults_",i,"_vs_dmea_absNES.pdf"),width=5,height=5)
  
  stats_spearman <- substitute(
    rho == est * "," ~ ~"p" ~ "=" ~ p,
    list(
      est = format(as.numeric(musyc.drug.corrsp$estimate), digits = 2),
      p = format(musyc.drug.corrsp$p.value, digits = 2)
    )
  )
  ggplot(shared.moa.res, aes(x=rank, y=abs(NES))) + # could set shape=sample but only sample is JH-2-002 
    geom_point(aes(color=moaCombo), show.legend =TRUE) + theme_minimal() + scale_color_manual(values=combo.cols, breaks=names(combo.cols),labels=names(combo.cols)) + 
    geom_smooth(method="lm",linetype="dashed",color="black") + ggplot2::geom_text(
      x = -Inf, y = Inf, vjust = "inward", hjust = "inward",
      colour = "black", parse = TRUE, 
      label = as.character(as.expression(stats_spearman)), size = 8
    ) +
    labs(x=i,
         y="Absolute DMEA Sensitivity Prediction")
  ggsave(paste0("moaResults_",i,"_vs_dmea_absNES_spearman.pdf"),width=5,height=5)
}

syn.scores <- c("Maximum Bliss Synergy Score" = "maxBliss", "Mean MuSyC Potency Score" = "meanLogAlpha")
shared.moa.res2 <- shared.moa.res[shared.moa.res$sig,]
for (i in names(syn.scores)) {
  shared.moa.res2$rank <- shared.moa.res2[,syn.scores[i]]
  # signed
  musyc.drug.corr <- cor.test(shared.moa.res2$rank, shared.moa.res2$NES)
  musyc.drug.corr$estimate # 0.2382369 
  musyc.drug.corr$p.value # 1.293807e-12
  musyc.drug.corrsp <- cor.test(shared.moa.res2$rank, shared.moa.res2$NES,method="spearman")
  musyc.drug.corrsp$estimate # 0.2414861
  musyc.drug.corrsp$p.value # 6.260826e-13
  
  stats_pearson <- substitute(
    r == est * "," ~ ~"p" ~ "=" ~ p,
    list(
      est = format(as.numeric(musyc.drug.corr$estimate), digits = 2),
      p = format(musyc.drug.corr$p.value, digits = 2)
    )
  )
  ggplot(shared.moa.res2, aes(x=rank, y=NES)) + # could set shape=sample but only sample is JH-2-002 
    geom_point(aes(color=moaCombo), show.legend =TRUE) + theme_minimal() + scale_color_manual(values=combo.cols, breaks=names(combo.cols),labels=names(combo.cols)) + 
    geom_smooth(method="lm",linetype="dashed",color="black") + ggplot2::geom_text(
      x = -Inf, y = Inf, vjust = "inward", hjust = "inward",
      colour = "black", parse = TRUE, 
      label = as.character(as.expression(stats_pearson)), size = 8
    ) +
    labs(x=i,
         y="DMEA Sensitivity Prediction")
  ggsave(paste0("moaResults_",i,"_vs_dmea_sigNES.pdf"),width=5,height=5)
  
  stats_spearman <- substitute(
    rho == est * "," ~ ~"p" ~ "=" ~ p,
    list(
      est = format(as.numeric(musyc.drug.corrsp$estimate), digits = 2),
      p = format(musyc.drug.corrsp$p.value, digits = 2)
    )
  )
  ggplot(shared.moa.res2, aes(x=rank, y=NES)) + # could set shape=sample but only sample is JH-2-002 
    geom_point(aes(color=moaCombo), show.legend =TRUE) + theme_minimal() + scale_color_manual(values=combo.cols, breaks=names(combo.cols),labels=names(combo.cols)) + 
    geom_smooth(method="lm",linetype="dashed",color="black") + ggplot2::geom_text(
      x = -Inf, y = Inf, vjust = "inward", hjust = "inward",
      colour = "black", parse = TRUE, 
      label = as.character(as.expression(stats_spearman)), size = 8
    ) +
    labs(x=i,
         y="DMEA Sensitivity Prediction")
  ggsave(paste0("moaResults_",i,"_vs_dmea_sigNES_spearman.pdf"),width=5,height=5)
  
  # abs drug corr
  musyc.drug.corr <- cor.test(shared.moa.res2$rank, abs(shared.moa.res2$NES))
  musyc.drug.corr$estimate # 0.2382369 
  musyc.drug.corr$p.value # 1.293807e-12
  musyc.drug.corrsp <- cor.test(shared.moa.res2$rank, abs(shared.moa.res2$NES),method="spearman")
  musyc.drug.corrsp$estimate # 0.2414861
  musyc.drug.corrsp$p.value # 6.260826e-13
  
  stats_pearson <- substitute(
    r == est * "," ~ ~"p" ~ "=" ~ p,
    list(
      est = format(as.numeric(musyc.drug.corr$estimate), digits = 2),
      p = format(musyc.drug.corr$p.value, digits = 2)
    )
  )
  ggplot(shared.moa.res2, aes(x=rank, y=abs(NES))) + # could set shape=sample but only sample is JH-2-002 
    geom_point(aes(color=moaCombo), show.legend =TRUE) + theme_minimal() + scale_color_manual(values=combo.cols, breaks=names(combo.cols),labels=names(combo.cols)) + 
    geom_smooth(method="lm",linetype="dashed",color="black") + ggplot2::geom_text(
      x = -Inf, y = Inf, vjust = "inward", hjust = "inward",
      colour = "black", parse = TRUE, 
      label = as.character(as.expression(stats_pearson)), size = 8
    ) +
    labs(x=i,
         y="Absolute DMEA Sensitivity Prediction")
  ggsave(paste0("moaResults_",i,"_vs_dmea_absSigNES.pdf"),width=5,height=5)
  
  stats_spearman <- substitute(
    rho == est * "," ~ ~"p" ~ "=" ~ p,
    list(
      est = format(as.numeric(musyc.drug.corrsp$estimate), digits = 2),
      p = format(musyc.drug.corrsp$p.value, digits = 2)
    )
  )
  ggplot(shared.moa.res2, aes(x=rank, y=abs(NES))) + # could set shape=sample but only sample is JH-2-002 
    geom_point(aes(color=moaCombo), show.legend =TRUE) + theme_minimal() + scale_color_manual(values=combo.cols, breaks=names(combo.cols),labels=names(combo.cols)) + 
    geom_smooth(method="lm",linetype="dashed",color="black") + ggplot2::geom_text(
      x = -Inf, y = Inf, vjust = "inward", hjust = "inward",
      colour = "black", parse = TRUE, 
      label = as.character(as.expression(stats_spearman)), size = 8
    ) +
    labs(x=i,
         y="Absolute DMEA Sensitivity Prediction")
  ggsave(paste0("moaResults_",i,"_vs_dmea_absSigNES_spearman.pdf"),width=5,height=5)
}

# top is MEK+HDAC 8h in MN-2 for DMEA, MEK+CDK for MuSyC


# ggplot(shared.moa.res, aes(x=meanLogAlpha, y=NES, #shape=sample, # only sample is JH-2-002
#                                                         color=moaComboShort)) +
#   geom_point() + theme_minimal() + #scale_color_manual(values=combo.cols, breaks=names(combo.cols), labels=names(combo.cols)) +
#   labs(x="MuSyC: % Change in Viability with Combination vs. Most Effective Single Agent",
#        y="DMEA: Predicted Drug Sensitivity")
# musyc.drug.corr <- cor.test(shared.moa.res$meanLogAlpha, shared.moa.res$NES)
# musyc.drug.corr$estimate # -0.07149301 
# musyc.drug.corr$p.value # 0.08647383
# musyc.drug.corrsp <- cor.test(shared.moa.res$meanLogAlpha, shared.moa.res$NES ,method="spearman")
# musyc.drug.corrsp$estimate # -0.09187795 
# musyc.drug.corrsp$p.value # 0.02745786
# 
# ggplot(shared.moa.res[shared.moa.res$sig,], aes(x=meanLogAlpha, y=NES, #shape=sample, # only sample is JH-2-002
#                            color=moaComboShort)) +
#   geom_point() + theme_minimal() + #scale_color_manual(values=combo.cols, breaks=names(combo.cols), labels=names(combo.cols)) +
#   labs(x="MuSyC: % Change in Viability with Combination vs. Most Effective Single Agent",
#        y="DMEA: Predicted Drug Sensitivity")
# musyc.drug.corr <- cor.test(shared.moa.res[shared.moa.res$sig,]$meanLogAlpha, shared.moa.res[shared.moa.res$sig,]$NES)
# musyc.drug.corr$estimate # 0.09917356
# musyc.drug.corr$p.value # 0.09298546
# musyc.drug.corrsp <- cor.test(shared.moa.res[shared.moa.res$sig,]$meanLogAlpha, shared.moa.res[shared.moa.res$sig,]$NES ,method="spearman")
# musyc.drug.corrsp$estimate # 0
# musyc.drug.corrsp$p.value # 1
# 
# ggplot(shared.moa.res2, aes(y=meanLogAlpha, x=NES)) + 
#   geom_point(aes(shape=moaComboShort, color=sig), size=3) + 
#   theme_minimal() + ggrepel::geom_text_repel(data=subset(shared.moa.res2,sig),
#                                              aes(label=paste(as.character(Time), "RNA,\n",as.character(timeD),"viability"))) + 
#   labs(y="MuSyC Synergy Score", shape = "Drug Combination", 
#        x="DMEA Sensitivity Score", color = "Significant\nDMEA Result") + 
#   scale_color_manual(values=c("red","grey"), breaks=c(TRUE, FALSE))
# ggsave("MOA_results_musyc_vs_dmea.pdf", width=7, height=4)
# 
# ggplot(shared.moa.res2, aes(y=meanLogAlpha, x=NES)) + 
#   geom_point(aes(shape=moaComboShort, color=sig), size=3) + 
#   theme_minimal() + 
#   labs(y="MuSyC Synergy Score", shape = "Drug Combination", 
#        x="DMEA Sensitivity Score", color = "Significant\nDMEA Result") + 
#   scale_color_manual(values=c("red","grey"), breaks=c(TRUE, FALSE))
# ggsave("MOA_results_musyc_vs_dmea_noLabels.pdf", width=7, height=4)
# 
# ggplot(shared.moa.res2, aes(y=meanLogAlpha, x=NES)) + 
#   geom_point(aes(shape=moaComboShort, color=sig), size=3) + 
#   geom_smooth(method="lm",se=FALSE,formula = y ~ x + I(x^2), color="darkgrey", linetype="dashed") +
#   theme_minimal() + 
#   labs(y="MuSyC Synergy Score", shape = "Drug Combination", 
#        x="DMEA Sensitivity Score", color = "Significant\nDMEA Result") +
#   scale_color_manual(values=c("red","grey"), breaks=c(TRUE, FALSE))
# ggsave("MOA_results_musyc_vs_dmea_noLabels_wQuadraticFit.pdf", width=7, height=4)
# 
# ggplot(shared.moa.res2, aes(y=meanLogAlpha, x=NES)) + 
#   geom_point(aes(shape=moaComboShort, color=sig), size=3) + 
#   geom_smooth(method="lm",se=FALSE,formula = y ~ x + I(x^2), color="darkgrey", linetype="dashed") +
#   theme_minimal() + ggrepel::geom_text_repel(data=subset(shared.moa.res2,sig),
#                                              aes(label=paste(as.character(Time), "RNA,\n",as.character(timeD),"viability"))) + 
#   labs(y="MuSyC Synergy Score", shape = "Drug Combination", 
#        x="DMEA Sensitivity Score", color = "Significant\nDMEA Result") +
#   scale_color_manual(values=c("red","grey"), breaks=c(TRUE, FALSE))
# ggsave("MOA_results_musyc_vs_dmea_wQuadraticFit.pdf", width=7, height=4)
# write.csv(shared.moa.res2,"MOA_results_musyc_vs_dmea.csv", row.names=FALSE)
# 
# ggplot(shared.moa.res2[shared.moa.res2$meanLogAlpha>0,], aes(y=meanLogAlpha, x=NES)) + 
#   geom_point(aes(shape=moaComboShort, color=sig), size=3) + 
#   theme_minimal() + ggrepel::geom_text_repel(data=subset(shared.moa.res2[shared.moa.res2$meanLogAlpha>0,],sig),
#                                              aes(label=paste(as.character(Time), "RNA,\n",as.character(timeD),"viability"))) + 
#   labs(y="MuSyC Synergy Score", shape = "Drug Combination", 
#        x="DMEA Sensitivity Score", color = "Significant\nDMEA Result") + 
#   scale_color_manual(values=c("red","grey"), breaks=c(TRUE, FALSE))
# ggsave("MOA_results_musyc_vs_dmea_positiveAlpha.pdf", width=7, height=4)
# 
# ggplot(shared.moa.res2[shared.moa.res2$meanLogAlpha>0,], aes(y=meanLogAlpha, x=NES)) + 
#   geom_point(aes(shape=moaComboShort, color=sig), size=3) + 
#   geom_smooth(method="lm",se=FALSE,formula = y ~ x + I(x^2), color="darkgrey", linetype="dashed") +
#   theme_minimal() + ggrepel::geom_text_repel(data=subset(shared.moa.res2[shared.moa.res2$meanLogAlpha>0,],sig),
#                                              aes(label=paste(as.character(Time), "RNA,\n",as.character(timeD),"viability"))) + 
#   labs(y="MuSyC Synergy Score", shape = "Drug Combination", 
#        x="DMEA Sensitivity Score", color = "Significant\nDMEA Result") +
#   scale_color_manual(values=c("red","grey"), breaks=c(TRUE, FALSE))
# ggsave("MOA_results_musyc_vs_dmea_wQuadraticFit_positiveAlpha.pdf", width=7, height=4)
# 
# ggplot(shared.moa.res2[shared.moa.res2$meanLogAlpha>0,], aes(y=meanLogAlpha, x=NES)) + 
#   geom_point(aes(shape=moaComboShort, color=sig), size=3) + 
#   geom_smooth(method="lm",se=FALSE,formula = y ~ x + I(x^2), color="darkgrey", linetype="dashed") +
#   theme_minimal() + 
#   labs(y="MuSyC Synergy Score", shape = "Drug Combination", 
#        x="DMEA Sensitivity Score", color = "Significant\nDMEA Result") +
#   scale_color_manual(values=c("red","grey"), breaks=c(TRUE, FALSE))
# ggsave("MOA_results_musyc_vs_dmea_wQuadraticFit_positiveAlpha_noLabels.pdf", width=7, height=4)
# 
# # are meanLogAlpha values higher in significant DMEA results?
# b.test <- t.test(shared.moa.res2[shared.moa.res2$sig,]$meanLogAlpha, shared.moa.res2[!shared.moa.res2$sig,]$meanLogAlpha, alternative="greater")
# b.test$p.value # 0.03112967
# b.test$estimate # mean of x mean of y 
# # 2.8980339 0.4124156 
# length(shared.moa.res2[shared.moa.res2$sig,]$meanLogAlpha) # 12
# length(shared.moa.res2[!shared.moa.res2$sig,]$meanLogAlpha) # 12
# 
# b.test <- t.test(shared.moa.res2[shared.moa.res2$sig & shared.moa.res2$meanLogAlpha>0,]$meanLogAlpha, 
#                  shared.moa.res2[!shared.moa.res2$sig & shared.moa.res2$meanLogAlpha>0,]$meanLogAlpha, alternative="greater")
# b.test$p.value # 0.06034868
# b.test$estimate # mean of x mean of y 
# # 4.179305  3.074639 
# length(shared.moa.res2[shared.moa.res2$sig & abs(shared.moa.res2$meanLogAlpha)<2,]$meanLogAlpha) # 3
# length(shared.moa.res2[!shared.moa.res2$sig & abs(shared.moa.res2$meanLogAlpha)<2,]$meanLogAlpha) # 1

#### bliss heatmap ####
maxAbsBliss <- max(abs(max.bliss.tested$maxBliss)) # 65.87242
minMaxBliss <- min(max.bliss.tested$maxBliss) # 0
ggplot(max.bliss.tested, aes(x=sample, y=reorder(drugCombo, maxBliss), fill=maxBliss)) + geom_tile(stat="identity") + 
  facet_wrap(.~paste0(time,"d"))+theme_classic() + scale_fill_gradient(low="grey",high="red",limits=c(minMaxBliss,maxAbsBliss)) + 
  theme(axis.title = element_blank(), axis.text.x = element_text(angle=45, vjust=1, hjust=1)) + labs(fill="Max Bliss")
ggsave("bliss_maxSynergyTested_heatmap.pdf", width=4.5,height=4)

minNonZeroBliss <- min(max.bliss.tested[max.bliss.tested$maxBliss != 0,]$maxBliss) # 0
ggplot(max.bliss.tested[max.bliss.tested$maxBliss != 0,], aes(x=sample, y=reorder(drugCombo, maxBliss), fill=log10(maxBliss))) + geom_tile(stat="identity") + 
  facet_wrap(.~paste0(time,"d"))+theme_classic() + scale_fill_gradient2(low="blue",mid="grey",high="red",limits=c(log10(minNonZeroBliss),log10(maxAbsBliss))) + 
  theme(axis.title = element_blank(), axis.text.x = element_text(angle=45, vjust=1, hjust=1)) + labs(fill="Log|Max Bliss|")
ggsave("bliss_maxSynergyTested_log10_heatmap.pdf", width=4.5,height=4)

max.bliss$maxBliss <- as.numeric(max.bliss$maxBliss)
max.bliss <- max.bliss[!is.na(max.bliss$maxBliss),]
maxAbsBliss <- max(abs(max.bliss$maxBliss)) # 65.87242
minMaxBliss <- min(max.bliss$maxBliss) # 0
ggplot(max.bliss, aes(x=sample, y=reorder(drugCombo, maxBliss), fill=maxBliss)) + geom_tile(stat="identity") + 
  facet_wrap(.~paste0(time,"d"))+theme_classic() + scale_fill_gradient(low="grey",high="red",limits=c(minMaxBliss,maxAbsBliss)) + 
  theme(axis.title = element_blank(), axis.text.x = element_text(angle=45, vjust=1, hjust=1)) + labs(fill="Max Bliss")
ggsave("bliss_maxSynergy_heatmap.pdf", width=4.5,height=4)

minNonZeroBliss <- min(max.bliss[max.bliss$maxBliss != 0,]$maxBliss) # 0
ggplot(max.bliss[max.bliss$maxBliss != 0,], aes(x=sample, y=reorder(drugCombo, maxBliss), fill=log10(maxBliss))) + geom_tile(stat="identity") + 
  facet_wrap(.~paste0(time,"d"))+theme_classic() + scale_fill_gradient2(low="blue",mid="grey",high="red",limits=c(log10(minNonZeroBliss),log10(maxAbsBliss))) + 
  theme(axis.title = element_blank(), axis.text.x = element_text(angle=45, vjust=1, hjust=1)) + labs(fill="Log|Max Bliss|")
ggsave("bliss_maxSynergy_log10_heatmap.pdf", width=4.5,height=4)

#### does synergy score correlate with -log(p_viability) ####
# load musyc scores
setwd(dataPath)
musyc.scores <- read.csv("Deconvolved_musyc_20250812.csv")
musyc.scores$drugCombo <- paste0(musyc.scores$drug1,"+",musyc.scores$drug2)
musyc.scores$mpnst <- musyc.scores$sample
musyc.scores$sample <- paste0(musyc.scores$mpnst,'_',musyc.scores$time*24,"hours")
musyc.scores$meanLogAlpha <- NA
for (i in 1:nrow(musyc.scores)) {
  musyc.scores$meanLogAlpha[i] <- mean(c(musyc.scores$log_alpha12[i], musyc.scores$log_alpha21[i])) 
}
musyc.scores$shortCombo <- paste0(tolower(substr(musyc.scores$drug1, 1,4)), "+",
                                  tolower(substr(musyc.scores$drug2, 1,4)))

# load bliss scores
setwd(paste0("curves_",Sys.Date()))
max.bliss.tested <- read.csv("maxBlissTested.csv")
max.bliss.tested$mpnst <- max.bliss.tested$sample
max.bliss.tested$sample <- paste0(max.bliss.tested$mpnst,'_',max.bliss.tested$time*24,"hours")
max.bliss.tested$drug1 <- sub("[+].*","",max.bliss.tested$drugCombo)
max.bliss.tested$drug2 <- sub(".*[+]","",max.bliss.tested$drugCombo)
max.bliss.tested$shortCombo <- paste0(tolower(substr(max.bliss.tested$drug1, 1,4)), "+",
                                  tolower(substr(max.bliss.tested$drug2, 1,4)))

# load p-values for single vs. combo
p.df <- read.csv("tTests_viability_combo_lessThan_single.csv")

# merge
p.musyc <- merge(p.df, musyc.scores, by=c("drugCombo","sample")) # 116 rows
p.bliss <- merge(p.df, max.bliss.tested, by=c("drugCombo","sample")) # 111 rows

# correlations: bliss
pb.corr <- cor.test(-log10(p.bliss$p),p.bliss$maxBliss) # r=0.1866132, p=0.02177
pb.corr
pb.corr.sp <- cor.test(-log10(p.bliss$p),p.bliss$maxBliss, method="spearman") # r=0.2005121, p=0.01356
pb.corr.sp
ggplot(p.bliss,aes(x=-log10(p),y=maxBliss))+geom_point(aes(#shape=sample,
                                                           color=drugCombo))+ 
  theme_minimal()+ geom_smooth(method="lm", linetype="dashed", colour="black")+
  ggrepel::geom_text_repel(data=subset(p.bliss, p<=0.05), aes(label=shortCombo), box.padding = 0.5) + 
  #scale_color_manual(values=cols,breaks=tolower(names(drug.info)))+
  labs(x="-Log(P-value)",y="Maximum Bliss Synergy",
       shape="MPNST") + geom_text(aes(x=4, y=60), label="r=0.37, p=6.6E-5", size=4)
ggsave("viability_-log10P_vs_maxBlissTested_noShape_sigComboLabel.pdf",width=7,height=5)

stats_spearman <- substitute(
  rho == est * "," ~ ~"p" ~ "=" ~ p,
  list(
    est = format(as.numeric(pb.corr.sp$estimate), digits = 2),
    p = format(pb.corr.sp$p.value, digits = 2)
  )
)
ggplot(p.bliss,aes(x=-log10(p),y=maxBliss))+geom_point(aes(#shape=sample,
  color=drugCombo)
                                                           )+ 
  theme_minimal()+ geom_smooth(method="lm", linetype="dashed", colour="black")+
  ggrepel::geom_text_repel(data=subset(p.bliss, p<=0.05), aes(label=shortCombo), box.padding = 0.5) + 
  #scale_color_manual(values=cols,breaks=tolower(names(drug.info)))+
  labs(x="-Log(P-value)",y="Maximum Bliss Synergy",
       shape="MPNST") + geom_text(aes(x=4, y=60), parse=TRUE, label=as.character(as.expression(stats_spearman)), size=4)
ggsave("viability_-log10P_vs_maxBlissTested_spearman_noShape_sigComboLabel.pdf",width=7,height=5)
write.csv(p.bliss,"p_Bliss.csv", row.names=FALSE)

# correlations: musyc
pm.corr <- cor.test(-log10(p.musyc$p),p.musyc$meanLogAlpha) # r=0.1250346, p=0.1199
pm.corr
pm.corr.sp <- cor.test(-log10(p.musyc$p),p.musyc$meanLogAlpha, method="spearman") # r=0.1607403, p = 0.0451
pm.corr.sp
ggplot(p.musyc,aes(x=-log10(p),y=meanLogAlpha))+geom_point(aes(#shape=sample,
  color=drugCombo))+ 
  theme_minimal()+ geom_smooth(method="lm", linetype="dashed", colour="black")+
  ggrepel::geom_text_repel(data=subset(p.musyc, p<=0.05), aes(label=shortCombo), box.padding = 0.5) + 
  #scale_color_manual(values=cols,breaks=tolower(names(drug.info)))+
  labs(x="-Log(P-value)",y="Mean MuSyC Potency",
       shape="MPNST") + geom_text(aes(x=4, y=10), label="r=0.195, p=3.57E-2", size=4)
ggsave("viability_-log10P_vs_musycMeanLogAlpha_noShape_sigComboLabel.pdf",width=7,height=5)

stats_spearman <- substitute(
  rho == est * "," ~ ~"p" ~ "=" ~ p,
  list(
    est = format(as.numeric(pm.corr.sp$estimate), digits = 2),
    p = format(pm.corr.sp$p.value, digits = 2)
  )
)
ggplot(p.musyc,aes(x=-log10(p),y=meanLogAlpha))+geom_point(aes(#shape=sample,
  color=drugCombo)
)+ 
  theme_minimal()+ geom_smooth(method="lm", linetype="dashed", colour="black")+
  #scale_color_manual(values=cols,breaks=tolower(names(drug.info)))+
  ggrepel::geom_text_repel(data=subset(p.musyc, p<=0.05), aes(label=shortCombo), box.padding = 0.5) + 
  labs(x="-Log(P-value)",y="Mean MuSyC Potency",
       shape="MPNST") + geom_text(aes(x=4, y=10), parse=TRUE, label=as.character(as.expression(stats_spearman)), size=4)
ggsave("viability_-log10P_vs_musycMeanLogAlpha_spearman_noShape_sigComboLabel.pdf",width=7,height=5)
write.csv(p.musyc, "p_musyc.csv", row.names=FALSE)

# musyc is better correlated
# what are the top musyc results
top.musyc <- musyc.scores[musyc.scores$meanLogAlpha > 0 & 
                            musyc.scores$log_alpha12 > 0 & musyc.scores$log_alpha21 > 0,
                          c("drugCombo","sample","meanLogAlpha","log_alpha12","log_alpha21")] # 130 rows
length(unique(top.musyc$drugCombo)) # 18 unique drug combos out of 27 tested
top.musyc$absDeltaMeanLogAlpha <- abs(top.musyc$log_alpha12 - top.musyc$log_alpha21)
top.musyc$fracAbsDelta <- top.musyc$absDeltaMeanLogAlpha / top.musyc$meanLogAlpha
best.musyc <- top.musyc[top.musyc$fracAbsDelta <= 0.25,] # 9 rows
length(unique(best.musyc$drugCombo)) # 8 unique drug combos out of 27 tested
write.csv(top.musyc,"positive_musyc.csv", row.names=FALSE)
write.csv(best.musyc,"positive_musyc_fracMax0.25.csv")

#### EGFRi for MEKi ####
drug.info <- read.csv("https://raw.githubusercontent.com/BelindaBGarana/DMEA/refs/heads/shiny-app/Inputs/PRISM_secondary-screen-replicate-treatment-info.csv")
red.drug.info <- dplyr::distinct(drug.info[,c("name","moa","target")])
targets <- na.omit(unique(red.drug.info$target))
targets <- unique(unlist(strsplit(targets, ", "))) # 1026

# get positive TF networks for when EGFRi was negatively enriched for MEKi
egfr.sens.mek <- moa.df[moa.df$sig & moa.df$NES < 0 & moa.df$Drug_set == "EGFR inhibitor" & moa.df$DrugTreatment %in% c("Mirdametinib","Selumetinib","Trametinib"),]
pos.TF.centr <- data.frame()
neg.TF.centr <- data.frame()
tf.path <- "/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/MPNST-PDX-MT/RNAseq/synapse_based_code/TF_networks"
for (i in 1:nrow(egfr.sens.mek)) {
  pos.filename <- file.path(tf.path,paste0(egfr.sens.mek$MPNST[i],"_",egfr.sens.mek$DrugTreatment[i],"_",egfr.sens.mek$Time[i],"_positive_TF_centrality.csv"))
  neg.filename <- file.path(tf.path,paste0(egfr.sens.mek$MPNST[i],"_",egfr.sens.mek$DrugTreatment[i],"_",egfr.sens.mek$Time[i],"_negative_TF_centrality.csv"))
  if (file.exists(pos.filename)) {
    pos.TF.centr <- rbind(pos.TF.centr, read.csv(pos.filename)) 
  }
  if (file.exists(neg.filename)) {
    neg.TF.centr <- rbind(neg.TF.centr, read.csv(neg.filename)) 
  }
}
write.csv(pos.TF.centr, "positive_MEKiTreated_TFnetworkCentrality_EGFRiSensitive.csv", row.names=FALSE)
write.csv(neg.TF.centr, "negative_MEKiTreated_TFnetworkCentrality_EGFRiSensitive.csv", row.names=FALSE)

pos.target.centr <- pos.TF.centr[tolower(pos.TF.centr$name) %in% tolower(targets),] # 366 rows
pos.target.centr$direction <- "positive"
neg.target.centr <- neg.TF.centr[tolower(neg.TF.centr$name) %in% tolower(targets),] # 526 rows
neg.target.centr$direction <- "negative"
target.centr <- rbind(pos.target.centr, neg.target.centr)
write.csv(target.centr, "drugTargets_MEKiTreated_TFnetworkCentrality_EGFRiSensitive.csv", row.names=FALSE)
target.mean.centr <- plyr::ddply(target.centr,.(name),summarize,
                                     meanPosCentrality = mean(eigen_centrality[direction=="positive"]),
                                 meanNegCentrality = mean(eigen_centrality[direction=="negative"]),
                                 meanCentrality = mean(eigen_centrality)) # 436
write.csv(target.mean.centr,"drugTargets_MEKiTreated_meanTFnetworkCentrality_EGFRiSensitive.csv", row.names=FALSE)
target.mean.centr.wInfo <- merge(target.mean.centr, red.drug.info, by.x="name",by.y="target") # 333 - only looks at drugs with 1 target
write.csv(target.mean.centr.wInfo,"drugTargetsSoloInfo_MEKiTreated_meanTFnetworkCentrality_EGFRiSensitive.csv", row.names=FALSE)
target.mean.centr.solo <- dplyr::distinct(target.mean.centr.wInfo[,colnames(target.mean.centr)]) # 145 individual targets
write.csv(target.mean.centr.solo,"drugTargetsSolo_MEKiTreated_meanTFnetworkCentrality_EGFRiSensitive.csv", row.names=FALSE)

#### aggregate results for DMEA ####
sens.moa <- moa.df[moa.df$sig & moa.df$NES < 0,]
mean.sens <- plyr::ddply(sens.moa, .(Drug_set), summarize,
                         meanNES = mean(NES),
                         medianNES = median(NES),
                         drugTreatments = paste0(unique(DrugTreatment), collapse=", "),
                         N_drugTreatments = length(unique(DrugTreatment)),
                         MPNSTs = paste0(unique(MPNST), collapse=", "),
                         N_MPNSTs = length(unique(MPNST)),
                         Times = paste0(unique(Time), collapse=", "),
                         N_Times = length(unique(Time)))
write.csv(mean.sens,"meanNegativeNES_DMEA.csv", row.names=FALSE)
ggplot(mean.sens,aes(x=medianNES, y=reorder(Drug_set,-medianNES))) + 
  geom_bar(stat="identity", fill="red") + 
  theme_classic() + scale_y_discrete(position="right") + 
  theme(axis.title.y=element_blank()) + 
  labs(x="Median NES", title="Drug MOAs Predicted to be Toxic")
ggsave("medianNegativeNES_DMEA_barPlot.pdf",width=4, height=7)

ggplot(mean.sens,aes(x=medianNES, y=N_drugTreatments)) + 
  geom_point() + ggrepel::geom_text_repel(aes(label=Drug_set)) + 
  theme_classic() + labs(x="Median NES", y="# of DrugTreatments", 
                         title="Drug MOAs Predicted to be Toxic")
ggsave("medianNegativeNES_DMEA_scatterPlot.pdf",width=4, height=4)

#### single agent curves ####
sing.conf <- read.table(file.path(dataPath,"mpnst_drug_response.tsv"), sep="\t", header=TRUE)
sing.conf$X <- FALSE
library(drc) # curve source: answer by greenjune: https://stackoverflow.com/questions/36780357/plotting-dose-response-curves-with-ggplot2-and-drc
# Sara shared these 2 links: https://stackoverflow.com/questions/68209998/plot-drc-using-ggplot; https://forum.posit.co/t/extract-out-points-in-dose-response-curve-produced-by-drc/159433
mean.sing <- plyr::ddply(sing.conf, .(DOSE, improve_sample_id, Drug, time), summarize,
                         meanGROWTH=mean(GROWTH),
                         sdGROWTH=sd(GROWTH))
mean.sing$timeD <- paste0(mean.sing$time/24,"d")
conf.plot <- ggplot(mean.sing, aes(x=DOSE, y=100*meanGROWTH, color=improve_sample_id)) +
  facet_grid(timeD ~ Drug) +
  geom_smooth(linetype="dashed", se=FALSE, method=drc::drm, method.args=list(fct=L.4(), se=FALSE)) +
  scale_x_continuous(transform="log10") + geom_point() + 
  geom_errorbar(aes(ymin=100*(meanGROWTH-sdGROWTH), ymax=100*(meanGROWTH+sdGROWTH))) +
  #ggtitle(paste0(d1," + ",d2))+scale_color_manual(values=temp.cols2) +
  theme_classic(base_size=12) + labs(x="Concentration (uM)", y = "% Viability", color="MPNST") +
  theme(plot.title=element_text(hjust=0.5,face="bold"), axis.text.x=element_text(angle=45, hjust=1, vjust=1))
conf.plot # something is off- looks like MN2 isn't relative viability, same for some MN4
ggsave("relativeViability_singleAgents.pdf", conf.plot, width=18, height=3)

top.sing <- mean.sing[mean.sing$Drug %in% c("TNO155","mirdametinib","vorinostat"),]
top.sing$Drug <- factor(top.sing$Drug, levels=c("TNO155","mirdametinib","vorinostat"))
conf.plot <- ggplot(top.sing, 
                    aes(x=DOSE, y=100*meanGROWTH, color=improve_sample_id)) +
  facet_grid(Drug ~ timeD) +
  geom_smooth(linetype="dashed", se=FALSE, method=drc::drm, method.args=list(fct=L.4(), se=FALSE)) +
  scale_x_continuous(transform="log10") + geom_point() + 
  geom_errorbar(aes(ymin=100*(meanGROWTH-sdGROWTH), ymax=100*(meanGROWTH+sdGROWTH))) +
  #ggtitle(paste0(d1," + ",d2))+scale_color_manual(values=temp.cols2) +
  theme_classic(base_size=12) + labs(x="Concentration (uM)", y = "% Viability", color="MPNST") +
  theme(plot.title=element_text(hjust=0.5,face="bold"), axis.text.x=element_text(angle=45, hjust=1, vjust=1))
conf.plot # something is off- looks like MN2 isn't relative viability, same for some MN4
ggsave("relativeViability_mirdametinibVorinostatTNO155.pdf", conf.plot, width=6, height=8)