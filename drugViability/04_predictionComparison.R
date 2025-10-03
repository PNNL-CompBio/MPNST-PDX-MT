# compare synergy scores to drug sensitivity predictions
# author: Belinda B. Garana
library(synapser);library(ggplot2);library(RColorBrewer)
setwd("~/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/MPNST-PDX-MT/drugViability")
dataPath <- "~/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/MPNST-PDX-MT/drugViability"
synapser::synLogin()

#results <- read.csv("results_viabilityBlissMusyc.csv")
results <- read.csv(synapser::synGet("syn68900322")$path)
drugs <- unique(c(results$drug1, results$drug2))

#### add delta AUC ####
dAUC <- read.csv(synapser::synGet("syn69978141")$path)
colnames(dAUC)[2:3] <- c("time","sample")
dAUC$time <- as.numeric(sub("h","",dAUC$time))/24
#dAUC$drugCombo[!(dAUC$drugCombo %in% results$drugCombo)] # 0
#dAUC$sample[!(dAUC$sample %in% results$sample)] # 0
results <- merge(results, dAUC, by=c("drugCombo","sample","time"))
combos <- unique(results$drugCombo)

#### compare musyc alpha to DMEA ####
##### drug corr #####
drug.corr <- read.csv(synapser::synGet("syn69910848")$path) # was syn65672266 which now doesn't include mpnst factor
drug.corr[drug.corr$DrugTreatment == "Trabectidin",]$DrugTreatment <- "Trabectedin"

shared.drugs <- drugs[tolower(drugs) %in% tolower(c(drug.corr$Drug, drug.corr$DrugTreatment))] # 9: mirdametinib, olaparib, TNO155, capmatinib, doxorubicin, palbociclib, vorinostat, irinotecan, trabectedin
shared.combos <- tidyr::expand_grid(shared.drugs, shared.drugs)
shared.combos$drugCombo <- paste0(shared.combos$shared.drugs...1,"+",shared.combos$shared.drugs...2)
shared.combos$drugCombo2 <- paste0(shared.combos$shared.drugs...2,"+",shared.combos$shared.drugs...1)
shared.combos <- unique(c(shared.combos$drugCombo, shared.combos$drugCombo2))
shared.combos <- shared.combos[shared.combos %in% combos] # 22 / 26 tested combos

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
synapser::synStore(synapser::File("drug_dmea_vs_synergy.csv", parent="syn68258288"))

# create colors using color in-between drug colors
drug.info <- list("palbociclib" = "CDK inhibitor", "ribociclib" = "CDK inhibitor",
                  "trabectedin" = "Chemotherapy", "ifosfamide" = "DNA alkylating agent",
                  "decitabine" = "DNMT inhibitor", "vorinostat" = "HDAC inhibitor",
                  "mirdametinib" = "MEK inhibitor", "selumetinib" = "MEK inhibitor",
                  "trametinib" = "MEK inhibitor", "capmatinib" = "MET inhibitor",
                  "olaparib" = "PARP inhibitor", "RMC4630" = "SHP2 inhibitor",
                  "TNO155" = "SHP2 inhibitor", "doxorubicin" = "TOP inhibitor",
                  "irinotecan" = "TOP inhibitor")
cols=c(RColorBrewer::brewer.pal(8,'Dark2'),RColorBrewer::brewer.pal(15-8,'Set2'))
names(cols) <- names(drug.info)
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
combo.cols <- unique(combo.cols) # 18
drug.corr.res <- dplyr::distinct(drug.corr.res)

syn.scores <- c("Maximum Bliss Synergy Score" = "maxBliss", 
                "Mean MuSyC Potency Score" = "meanLogAlpha",
                "Delta AUC" = "dAUC")
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
moa.df <- read.csv(synapser::synGet("syn69910844")$path) # was syn65672258 which now doesn't include mpnst factor
if (any(moa.df$DrugTreatment=="Trabectidin")) {
  moa.df[moa.df$DrugTreatment=="Trabectidin",]$DrugTreatment <- "Trabectedin"
}
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
for (d in shared.drugs[shared.drugs %in% shared.moa.df$DrugTreatment]) {
  shared.moa.df[tolower(shared.moa.df$DrugTreatment) == tolower(d),]$moa <- drug.info[[d]]
}
length(unique(shared.moa.df$moa)) # 8
shared.moas2 <- shared.moas[shared.moas %in% shared.moa.df$moa] # 8

results2 <- results[,c("sample","drug1","drug2","meanLogAlpha","maxBliss","R2","time","dAUC")]
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
shared.moa.res <- dplyr::distinct(shared.moa.res) # 96 rows
write.csv(shared.moa.res, "moa_dmea_vs_synergy.csv", row.names = FALSE)
synapser::synStore(synapser::File("moa_dmea_vs_synergy.csv", parent="syn68258288"))

syn.scores <- c("Maximum Bliss Synergy Score" = "maxBliss", 
                "Mean MuSyC Potency Score" = "meanLogAlpha",
                "Delta AUC" = "dAUC")
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

syn.scores <- c("Maximum Bliss Synergy Score" = "maxBliss", 
                "Mean MuSyC Potency Score" = "meanLogAlpha",
                "Delta AUC" = "dAUC")
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
