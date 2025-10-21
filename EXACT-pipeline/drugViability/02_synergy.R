# compile synergy results
library(plyr);library(dplyr);library(ggplot2);library(ggrepel)
library(reshape2);library(RColorBrewer);library(ggcorrplot)
library(fmsb);library(stats);library(synapser)
# author: Belinda B. Garana
#setwd("~/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/MPNST-PDX-MT/drugViability")

synapser::synLogin()

#rel.conf <- read.csv("mpnst_combo_drug_response.csv")
rel.conf <- read.csv(synapser::synGet("syn68156852")$path)
musyc.scores <- read.csv(synapser::synGet("syn68736713")$path) # was: musyc_20250812.csv

# try looking at Bliss scores
bliss <- data.frame()

#dir.create("bliss")
##SG Alex uploaded a file we dont need to concatenate
bliss <- read.csv(synapser::synGet('syn68900210')$path)
# blissFiles <- synapser::as.list(synapser::synGetChildren('syn68639935', list("file"), sortBy = 'NAME'))
# for (i in 1:length(blissFiles)){
#   test.bliss <- read.csv(synapser::synGet(blissFiles[[i]]$id)$path)
#   bliss <- rbind(bliss, test.bliss)
# }
# bliss$drugCombo <- paste0(bliss$drug1, "+", bliss$drug2)
# bliss$sample <- bliss$PDX
# bliss$time <- bliss$timePoint..hr./24
# write.csv(bliss, "bliss.csv", row.names=FALSE)
#bliss <- read.csv("bliss.csv")

# calculate mean and sd for each drug & dose combo, mpnst, time combo; also preserve sample and chr8q columns
mean.conf <- plyr::ddply(rel.conf, .(sample, drug1, drug1.conc, drug2, drug2.conc), summarize,
                         meanGROWTH = mean(effect, na.rm=TRUE),
                         sdGROWTH = sd(effect, na.rm=TRUE))
mean.conf <- mean.conf |>
  tidyr::separate(sample, into = c('sample','time'), sep = '_')

#mean.conf$time <- sub(".*_","",mean.conf$sample)
mean.conf$time <- sub("hours","",mean.conf$time)
mean.conf$time <- as.numeric(mean.conf$time)/24

#musyc.scores$time <- sub(".*_","",musyc.scores$sample)
#musyc.scores$time <- sub("hours","",musyc.scores$time)
#musyc.scores$time <- as.numeric(musyc.scores$time)/24

# # fix sample names - no longer need to since no longer shortening names for musyc
# full.names <- c("MN-4","JH-2-079c", "JH-2-002", "WU-225","WU-356", "MN-2")
# short.names <- c("N-4","79c", "002", "225", "356", "N-2")
# for (i in 1:length(short.names)) {
#   mean.conf[grepl(short.names[[i]],mean.conf$sample),]$sample <- full.names[[i]]
#   musyc.scores[grepl(short.names[[i]],musyc.scores$sample),]$sample <- full.names[[i]]
# }
# prc2.loss <- c("JH-2-079c", "JH-2-002", "WU-225", "WU-365", "MN-2")
# prc2.wt <- c("MN-4")
# mean.conf$PRC2 <- "Deficient"
# mean.conf[mean.conf$sample %in% prc2.wt,]$PRC2 <- "Intact"
#
# # restore full drug names for musyc data
# drugs <- unique(c(mean.conf$drug1, mean.conf$drug2))
# short.drugs <- substr(drugs, 2, 4)
# colnames(musyc.scores)[2:3] <- c("drug1","drug2")
# for (i in 1:length(drugs)) {
#   if (any(musyc.scores$drug1 == short.drugs[[i]])) {
#     musyc.scores[musyc.scores$drug1 == short.drugs[[i]],]$drug1 <- drugs[[i]]
#   }
#   if (any(musyc.scores$drug2 == short.drugs[[i]])) {
#     musyc.scores[musyc.scores$drug2 == short.drugs[[i]],]$drug2 <- drugs[[i]]
#   }
# }
# write.csv(musyc.scores, "Deconvolved_musyc_20250812.csv", row.names=FALSE) # duplicated and changed name to musyc.csv for synapse
# #musyc.scores <- read.csv("Deconvolved_musyc_20250812.csv")
musyc.scores <- read.csv(synapser::synGet("syn68736713")$path) |>
  tidyr::separate(sample,into=c('sample','time'),sep='_')

musyc.scores$time <- sapply(musyc.scores$time,function(x) ifelse(x=='48hours',2,5))

# make scatter plot of a12 vs a21
dir.create(paste0("curves_",Sys.Date()))
setwd(paste0("curves_",Sys.Date()))
musyc.scores$drugCombo <- paste0(musyc.scores$drug1,"+",musyc.scores$drug2)
# ggplot(musyc.scores, aes(x=log_alpha21, y=log_alpha12, shape=sample, color=drugCombo, size=R2)) +
#   geom_point()+theme_minimal() +
#   #ggrepel::geom_text_repel(aes(label=paste0(time,"d"))) +
#   labs(x="Drug 2's Effect on Potency of Drug 1", y="Drug 2's Effect on Potency of Drug 1")
#
# ggplot(musyc.scores[musyc.scores$log_alpha12>0 & musyc.scores$log_alpha21>0,],
#        aes(x=log_alpha21, y=log_alpha12, shape=sample, color=drugCombo, size=R2)) +
#   geom_point()+theme_minimal() +
#   #ggrepel::geom_text_repel(aes(label=paste0(time,"d"))) +
#   labs(x="Drug 2's Effect on Potency of Drug 1", y="Drug 2's Effect on Potency of Drug 1")
#
# musyc.scores$shortCombo <- paste0(tolower(substr(musyc.scores$drug1, 1,4)), "+",
#                                   tolower(substr(musyc.scores$drug2, 1,4)))
# ggplot(musyc.scores[musyc.scores$log_alpha12>0 & musyc.scores$log_alpha21>0,],
#        aes(x=log_alpha21, y=log_alpha12, shape=as.factor(time), color=drugCombo, size=R2)) +
#   geom_point()+theme_classic() + facet_wrap(.~sample)+ scale_size_continuous(breaks=c(0.3,0.5,0.7)) +
#   ggrepel::geom_text_repel(aes(label=shortCombo), box.padding = 0.5) +
#   labs(x="Drug 2's Effect on Potency of Drug 1", y="Drug 2's Effect on Potency of Drug 1", shape="Time (d)")
# ggsave("musyc_potentCombos_v2.pdf",width=12, height=9)
# #ggsave("musyc_potentCombos_classicTheme.pdf",width=12, height=6)

musyc.scores$meanLogAlpha <- NA
for (i in 1:nrow(musyc.scores)) {
  musyc.scores$meanLogAlpha[i] <- mean(c(musyc.scores$log_alpha12[i], musyc.scores$log_alpha21[i]))
}
maxLogAlpha <- max(musyc.scores$meanLogAlpha)
musyc.scores$drugCombo <- paste0(musyc.scores$drug1_name,'+',musyc.scores$drug2_name)
##remove combos requested by alex (SG)
musyc.scores <- subset(musyc.scores,!drugCombo%in%c("selumetinib+nazartinib","mirdametinib+nazartinib","avutometinib+defactinib"))
# ggplot(musyc.scores, aes(x=sample, y=reorder(drugCombo, meanLogAlpha), fill=meanLogAlpha)) + geom_tile(stat="identity") +
#   facet_wrap(.~paste0(time,"d"))+theme_classic() + scale_fill_gradient(low="grey",high="red",limits=c(0,8)) +
#   theme(axis.title = element_blank(), axis.text.x = element_text(angle=45, vjust=1, hjust=1)) + labs(fill="Mean\nLog|Alpha|")
# ggsave("musyc_meanLogAlpha_heatmap_positive.pdf", width=4,height=4)

ggplot(musyc.scores, aes(x=sample, y=reorder(drugCombo, meanLogAlpha), fill=meanLogAlpha)) + geom_tile(stat="identity") +
  facet_wrap(.~paste0(time,'d'))+theme_classic() + scale_fill_gradient2(low="blue",mid="grey",high="red",limits=c(-8,8)) +
  theme(axis.title = element_blank(), axis.text.x = element_text(angle=45, vjust=1, hjust=1)) + labs(fill="Mean\nLog|Alpha|")
ggsave("musyc_meanLogAlpha_heatmap.pdf", width=5,height=4)

 bliss.musyc <- merge(bliss, musyc.scores, by=c("drugCombo","sample","time"))
 max.bliss <- plyr::ddply(bliss, .(drugCombo, sample, time), summarize,
                          maxBliss = max(Bliss_synergy),
                          medianBliss = median(Bliss_synergy)) # 125
# write.csv(max.bliss,"maxBliss.csv", row.names=FALSE)

mean.conf$conc1 <- mean.conf$drug1.conc
mean.conf$conc2 <- mean.conf$drug2.conc
bliss.tested <- merge(bliss, mean.conf, by=c("sample","time","drug1","drug2","conc1","conc2")) # 2529 down from 10287 bliss
max.bliss.tested <- plyr::ddply(bliss.tested, .(drugCombo, sample, time), summarize,
                                maxBliss = max(Bliss_synergy),
                                medianBliss = median(Bliss_synergy),
                                meanBliss = mean(Bliss_synergy)) # 111
write.csv(max.bliss.tested,"maxBlissTested.csv", row.names=FALSE)

# bliss.musyc <- merge(max.bliss, musyc.scores, by=c("drugCombo","sample","time"))
# ggplot(bliss.musyc, aes(x=meanLogAlpha, y=maxBliss, shape=sample, color=drugCombo)) + #, size=R2)) +
#   geom_point()+theme_minimal() +
#   #ggrepel::geom_text_repel(aes(label=paste0(time,"d"))) +
#   labs(x="Mean MuSyC Potency", y="Max Bliss Synergy")
# ggsave("musyc_meanLogAlpha_vs_max_bliss_synergy.pdf", width=10, height=8)
#
# ggplot(bliss.musyc, aes(x=meanLogAlpha, y=medianBliss, shape=sample, color=drugCombo)) + #, size=R2)) +
#   geom_point()+theme_minimal() +
#   #ggrepel::geom_text_repel(aes(label=paste0(time,"d"))) +
#   labs(x="Mean MuSyC Potency", y="Max Bliss Synergy")
# ggsave("musyc_meanLogAlpha_vs_median_bliss_synergy.pdf", width=10, height=8)
#
# ggplot(bliss.musyc, aes(x=meanLogAlpha, y=maxBliss, shape=sample, color=drugCombo)) + #, size=R2)) +
#   geom_point()+theme_minimal() + facet_grid(paste0(time,"d")~sample)+
#   #ggrepel::geom_text_repel(aes(label=paste0(time,"d"))) +
#   labs(x="Mean MuSyC Potency", y="Max Bliss Synergy")
# ggsave("musyc_meanLogAlpha_vs_max_bliss_synergy_Faceted.pdf", width=12, height=7)
#
# ggplot(bliss.musyc, aes(x=meanLogAlpha, y=medianBliss, shape=sample, color=drugCombo)) + #, size=R2)) +
#   geom_point()+theme_minimal() + facet_grid(paste0(time,"d")~sample)+
#   #ggrepel::geom_text_repel(aes(label=paste0(time,"d"))) +
#   labs(x="Mean MuSyC Potency", y="Max Bliss Synergy")
# ggsave("musyc_meanLogAlpha_vs_median_bliss_synergy_Faceted.pdf", width=12, height=7)
#
# cor.test(bliss.musyc$meanLogAlpha, bliss.musyc$maxBliss) # r = 0.011, p=0.89
# cor.test(bliss.musyc$meanLogAlpha, bliss.musyc$maxBliss, method="spearman") # rho = 0.039, p = 0.6151
#
# cor.test(bliss.musyc$meanLogAlpha, bliss.musyc$medianBliss) # r = 0.1604656, p=0.0389
# cor.test(bliss.musyc$meanLogAlpha, bliss.musyc$medianBliss, method="spearman") # rho = 0.109585, p = 0.1599

bliss.musyc.tested <- merge(max.bliss.tested, musyc.scores, by=c("drugCombo","sample","time"))
# ggplot(bliss.musyc.tested, aes(x=meanLogAlpha, y=maxBliss, shape=sample, color=drugCombo)) + #, size=R2)) +
#   geom_point()+theme_minimal() +
#   #ggrepel::geom_text_repel(aes(label=paste0(time,"d"))) +
#   labs(x="Mean MuSyC Potency", y="Max Bliss Synergy")
# ggsave("musyc_meanLogAlpha_vs_max_bliss_synergy_tested.pdf", width=10, height=8)
#
# ggplot(bliss.musyc.tested, aes(x=meanLogAlpha, y=maxBliss, shape=sample, color=drugCombo)) + #, size=R2)) +
#   geom_point()+theme_minimal() + facet_grid(paste0(time,"d")~sample)+
#   #ggrepel::geom_text_repel(aes(label=paste0(time,"d"))) +
#   labs(x="Mean MuSyC Potency", y="Max Bliss Synergy")
# ggsave("musyc_meanLogAlpha_vs_max_bliss_synergy_tested_Faceted.pdf", width=12, height=7)
#
# cor.test(bliss.musyc.tested$meanLogAlpha, bliss.musyc.tested$maxBliss) # r=0.066, p=0.3994
# cor.test(bliss.musyc.tested$meanLogAlpha, bliss.musyc.tested$maxBliss, method="spearman") # rho = 0.1481161, p=0.0561
#
# ggplot(bliss.musyc.tested, aes(x=meanLogAlpha, y=medianBliss, shape=sample, color=drugCombo)) + #, size=R2)) +
#   geom_point()+theme_minimal() +
#   #ggrepel::geom_text_repel(aes(label=paste0(time,"d"))) +
#   labs(x="Mean MuSyC Potency", y="Median Bliss Synergy")
# ggsave("musyc_meanLogAlpha_vs_median_bliss_synergy_tested.pdf", width=10, height=8)
#
# ggplot(bliss.musyc.tested, aes(x=meanLogAlpha, y=medianBliss, shape=sample, color=drugCombo)) + #, size=R2)) +
#   geom_point()+theme_minimal() + facet_grid(paste0(time,"d")~sample)+
#   #ggrepel::geom_text_repel(aes(label=paste0(time,"d"))) +
#   labs(x="Mean MuSyC Potency", y="Median Bliss Synergy")
# ggsave("musyc_meanLogAlpha_vs_median_bliss_synergy_tested_Faceted.pdf", width=12, height=7)
#
# cor.test(bliss.musyc.tested$meanLogAlpha, bliss.musyc.tested$medianBliss) # r=0.066, p=0.3994
# cor.test(bliss.musyc.tested$meanLogAlpha, bliss.musyc.tested$medianBliss, method="spearman") # rho = 0.1481161, p=0.0561

#### bliss heatmap ####
maxAbsBliss <- max(abs(max.bliss.tested$maxBliss)) # 65.87242
minMaxBliss <- min(max.bliss.tested$maxBliss) # 0
ggplot(max.bliss.tested, aes(x=sample, y=reorder(drugCombo, maxBliss), fill=maxBliss)) + geom_tile(stat="identity") +
  facet_wrap(.~paste0(time,"d"))+theme_classic() + scale_fill_gradient(low="grey",high="red",limits=c(minMaxBliss,maxAbsBliss)) +
  theme(axis.title = element_blank(), axis.text.x = element_text(angle=45, vjust=1, hjust=1)) + labs(fill="Max Bliss")
ggsave("bliss_maxSynergyTested_heatmap.pdf", width=4.5,height=4)

# minNonZeroBliss <- min(max.bliss.tested[max.bliss.tested$maxBliss != 0,]$maxBliss) # 0
# ggplot(max.bliss.tested[max.bliss.tested$maxBliss != 0,], aes(x=sample, y=reorder(drugCombo, maxBliss), fill=log10(maxBliss))) + geom_tile(stat="identity") +
#   facet_wrap(.~paste0(time,"d"))+theme_classic() + scale_fill_gradient2(low="blue",mid="grey",high="red",limits=c(log10(minNonZeroBliss),log10(maxAbsBliss))) +
#   theme(axis.title = element_blank(), axis.text.x = element_text(angle=45, vjust=1, hjust=1)) + labs(fill="Log|Max Bliss|")
# ggsave("bliss_maxSynergyTested_log10_heatmap.pdf", width=4.5,height=4)
#

##NEXT CHUNK commented by SG
 max.bliss$maxBliss <- as.numeric(max.bliss$maxBliss)
 max.bliss <- max.bliss[!is.na(max.bliss$maxBliss),]
 maxAbsBliss <- max(abs(max.bliss$maxBliss)) # 65.87242
 minMaxBliss <- min(max.bliss$maxBliss) # 0
 ggplot(max.bliss, aes(x=sample, y=reorder(drugCombo, maxBliss), fill=maxBliss)) + geom_tile(stat="identity") +
   facet_wrap(.~paste0(time,"d"))+theme_classic() + scale_fill_gradient(low="grey",high="red",limits=c(minMaxBliss,maxAbsBliss)) +
   theme(axis.title = element_blank(), axis.text.x = element_text(angle=45, vjust=1, hjust=1)) + labs(fill="Max Bliss")
 ggsave("bliss_maxSynergy_heatmap.pdf", width=4.5,height=4)

# minNonZeroBliss <- min(max.bliss[max.bliss$maxBliss != 0,]$maxBliss) # 0
# ggplot(max.bliss[max.bliss$maxBliss != 0,], aes(x=sample, y=reorder(drugCombo, maxBliss), fill=log10(maxBliss))) + geom_tile(stat="identity") +
#   facet_wrap(.~paste0(time,"d"))+theme_classic() + scale_fill_gradient2(low="blue",mid="grey",high="red",limits=c(log10(minNonZeroBliss),log10(maxAbsBliss))) +
#   theme(axis.title = element_blank(), axis.text.x = element_text(angle=45, vjust=1, hjust=1)) + labs(fill="Log|Max Bliss|")
# ggsave("bliss_maxSynergy_log10_heatmap.pdf", width=4.5,height=4)

maxAbsBliss <- max(abs(max.bliss.tested$medianBliss)) # 0
minMaxBliss <- min(max.bliss.tested$medianBliss) # 0
ggplot(max.bliss.tested, aes(x=sample, y=reorder(drugCombo, medianBliss), fill=medianBliss)) + geom_tile(stat="identity") +
  facet_wrap(.~paste0(time,"d"))+theme_classic() + scale_fill_gradient(low="grey",high="red",limits=c(minMaxBliss,maxAbsBliss)) +
  theme(axis.title = element_blank(), axis.text.x = element_text(angle=45, vjust=1, hjust=1)) + labs(fill="Median Bliss")
ggsave("bliss_medianSynergyTested_heatmap.pdf", width=4.5,height=4)

maxAbsBliss <- max(abs(max.bliss.tested$meanBliss)) # 24.871
minMaxBliss <- min(max.bliss.tested$meanBliss) # -24.871
ggplot(max.bliss.tested, aes(x=sample, y=reorder(drugCombo, meanBliss), fill=meanBliss)) + geom_tile(stat="identity") +
  facet_wrap(.~paste0(time,"d"))+theme_classic() + scale_fill_gradient2(low="blue",mid="grey",high="red",limits=c(minMaxBliss,maxAbsBliss)) +
  theme(axis.title = element_blank(), axis.text.x = element_text(angle=45, vjust=1, hjust=1)) + labs(fill="Mean Bliss")
ggsave("bliss_meanSynergyTested_heatmap.pdf", width=4.5,height=4)

#### spider plots ####
# 6 drug combos tested in all 6 MPNSTs at 2 timepoints: mean MuSyC LogAlpha, median CI, median Bliss
library(fmsb)

# get CI
ci.48h <- read.csv(synapser::synGet("syn70365485")$path) # 48hr Median CI of experimental doses.csv
ci.120h <- read.csv(synapser::synGet("syn70365484")$path) # 120hr Median CI of experimental doses.csv
ci.48h$time <- 2
ci.120h$time <- 5
ci <- rbind(ci.48h, ci.120h)
colnames(ci)[1] <- "drugCombo"
ci$drugCombo <- tolower(ci$drugCombo)
ci <- reshape2::melt(ci, id.vars=c("drugCombo","time"), variable.name="sample", value.name="CI")|>
  subset(!is.na(CI))
ci$sample <- gsub("[.]","-",ci$sample)

# cap abs(CI) values at 99 per Alex Larsson
#ci[ci$CI > 99,]$CI <- 99
#ci[ci$CI < -99,]$CI <- -99

# switch direction so that positive values correspond to synergy across all metrics
ci$rank <- -ci$CI # ci < 1 denotes synergy, =1: additive, <1: antagonistic
ci$CI <- NULL

# get delta AUC
#SG updated to synapse
dauc <- read.csv(syn$get('syn69978141')$path)##file.path(dataPath,"dAUC.csv"))
colnames(dauc) <- colnames(ci)
dauc$time <- as.numeric(sub("h","",dauc$time))/24

maxmin <- data.frame(minDauc = floor(min(dauc$rank, na.rm=TRUE)),maxDauc = ceiling(max(dauc$rank, na.rm=TRUE)),
                     minCI = floor(min(ci$rank, na.rm=TRUE)),maxCI = ceiling(max(ci$rank, na.rm=TRUE)),
                     minBliss = floor(min(bliss.musyc.tested$maxBliss, na.rm=TRUE)), maxBliss = ceiling(max(bliss.musyc.tested$maxBliss, na.rm=TRUE)),
                     minMusyc = floor(min(bliss.musyc.tested$meanLogAlpha, na.rm=TRUE)), maxMusyc = ceiling(max(bliss.musyc.tested$meanLogAlpha, na.rm=TRUE)))
maxmin2 <- data.frame("Delta AUC" = c(maxmin$maxDauc, maxmin$minDauc),
                      "CI" = c(maxmin$maxCI, maxmin$minCI),
                      "Bliss" = c(maxmin$maxBliss, maxmin$minBliss),
                      "MuSyC" = c(maxmin$maxMusyc, maxmin$minMusyc))
rownames(maxmin2) <- c("max","min")
bliss.musyc.tested$drugCombo <- tolower(bliss.musyc.tested$drugCombo)
ci$drugCombo <- tolower(ci$drugCombo)
synergy <- merge(na.omit(ci), bliss.musyc.tested[,c("drugCombo","sample","time","maxBliss","meanLogAlpha")],
                 by=c("drugCombo","sample","time"))
colnames(dauc)[ncol(dauc)] <- "dAUC"
synergy <- merge(na.omit(dauc), synergy,
                 by=c("drugCombo","sample","time"))
# comboCoverage <- plyr::ddply(synergy, .(drugCombo), summarize,
#                              nSamples=length(unique(sample))) # only mirda+palbo, mirda+vorin, mirda+trab have full coverage (6 samples); NAs at 5d CI in TNO combos
# maxSamples <- length(unique(synergy$sample)) # 6
# synergyFull <- synergy[synergy$drugCombo %in% unique(comboCoverage[comboCoverage$nSamples==maxSamples,]$drugCombo),]
# synergyFull <- reshape2::melt(synergyFull, id.vars=c("drugCombo","sample","time"), variable.name="Metric", value.name="Score")
# synergyFull <- reshape2::dcast(synergyFull, Metric+sample+time ~ drugCombo, value.var="Score") # there are NAs for MN-4 in mirda+palbo and mirda+vorin... make plots for combos instead of metrics

synergyFull <- reshape2::melt(synergy, id.vars=c("drugCombo","sample","time"), variable.name="Metric", value.name="Score")
synergyFull <- reshape2::dcast(synergyFull, drugCombo+sample+time ~ Metric, value.var="Score")
allSampleColors <- RColorBrewer::brewer.pal(length(unique(synergyFull$sample)),"Set2")
names(allSampleColors) <- unique(synergyFull$sample)

for (t in unique(synergyFull$time)) {
  for (dc in unique(synergyFull$drugCombo)) {
    # reformat data so that each column is an axis variable, rownames are samples except for first 2 rows which are max and min
    tempSynergy <- synergyFull[synergyFull$time == t & synergyFull$drugCombo == dc,]
    rownames(tempSynergy) <- tempSynergy$sample
    data <- tempSynergy[,c("dAUC","rank","maxBliss","meanLogAlpha")]
    colnames(data) <- c("Delta AUC","CI", "Bliss", "MuSyC")
    #data <- rbind(maxmin2, tempSynergy)
    sampleColors <- allSampleColors[rownames(data)]
    if (length(unique(data$Bliss))>1 & length(unique(data$CI)) > 1 & length(unique(data$MuSyC)) > 1) {
      # prep to save plot
      pdf(paste0(t,"d_",dc,"_spiderPlot_",Sys.Date(),".pdf"), width=5, height=4)

      # create plot
      fmsb::radarchart(data, axistype=0, maxmin=FALSE, pcol=sampleColors,
                       pfcol=alpha(sampleColors,0.3), plwd=4, plty=1, cglcol="gray",
                       cglty=1, axislabcol="black", cglwd=0.8, vlcex=0.8)
      legend(x=0.75, y=1.25, legend=rownames(data), bty="n", pch=20,
             col=sampleColors, text.col="black", cex=0.8, pt.cex=1.2)

      # save plot
      dev.off()
    }
  }
}

#### correlation matrices ####
# correlation matrix - Pearson
med.mat <- synergy
med.mat$sample <- paste0(med.mat$drugCombo,"_",med.mat$sample,"_",med.mat$time,"d")
rownames(med.mat) <- med.mat$sample
colnames(med.mat)[4:7] <- c("Delta AUC","CI","Bliss","MuSyC")
med.mat <- as.matrix(med.mat[,4:7])
corr.mat <- stats::cor(med.mat)
write.csv(corr.mat,"synergy_PearsonCorr.csv")

# plot correlation matrix
corr.mat.plot <- ggcorrplot::ggcorrplot(corr.mat)
ggsave("synergy_PearsonCorr.pdf", corr.mat.plot, width=4, height=4)

# correlation matrix - Spearman
corr.mat <- stats::cor(med.mat, method="spearman")
write.csv(corr.mat,"synergy_SpearmanCorr.csv")

# plot correlation matrix
corr.mat.plot <- ggcorrplot::ggcorrplot(corr.mat)
ggsave("synergy_SpearmanCorr.pdf", corr.mat.plot, width=4, height=4)

#### individual correlations ####
others <- c("CI","Bliss","MuSyC")
moa.med <- synergy
colnames(moa.med)[4:7] <- c("Delta AUC", others)
for (i in others) {
  moa.med$rank <- moa.med[,i]
  moa.med.cor <- cor.test(moa.med$`Delta AUC`, moa.med$rank)
  stats_pearson <- substitute(
    r == est * "," ~ ~"p" ~ "=" ~ p,
    list(
      est = format(as.numeric(moa.med.cor$estimate), digits = 3),
      p = format(moa.med.cor$p.value, digits = 3)
    )
  )
  ggplot(moa.med, aes(x=`Delta AUC`, y=rank)) + # could set shape=sample but only sample is JH-2-002
    geom_point(aes(color=drugCombo, shape=sample, size=as.factor(time))) + theme_minimal() +
    geom_smooth(method="lm",linetype="dashed",color="black") + ggplot2::geom_text(
      x = -Inf, y = Inf, vjust = "inward", hjust = "inward",
      colour = "black", parse = TRUE,
      label = as.character(as.expression(stats_pearson)), size = 6
    ) +
    labs(x="Delta AUC",color="Drug Combo", shape="MPNST", size="Time (days)",
         y=i)
  ggsave(paste0(i,"_vs_","dAUC","_Pearson.pdf"),width=7,height=5)

  moa.med.cor <- cor.test(moa.med$`Delta AUC`, moa.med$rank, method="spearman")
  stats_spearman <- substitute(
    rho == est * "," ~ ~"p" ~ "=" ~ p,
    list(
      est = format(as.numeric(moa.med.cor$estimate), digits = 3),
      p = format(moa.med.cor$p.value, digits = 3)
    )
  )
  ggplot(moa.med, aes(x=`Delta AUC`, y=rank)) + # could set shape=sample but only sample is JH-2-002
    geom_point(aes(color=drugCombo, shape=sample, size=as.factor(time))) + theme_minimal() +
    geom_smooth(method="lm",linetype="dashed",color="black") + ggplot2::geom_text(
      x = -Inf, y = Inf, vjust = "inward", hjust = "inward",
      colour = "black", parse = TRUE,
      label = as.character(as.expression(stats_spearman)), size = 6
    ) +
    labs(x="Delta AUC",color="Drug Combo", shape="MPNST", size="Time (days)",
         y=i)
  ggsave(paste0(i,"_vs_","dAUC","_Spearman.pdf"),width=7,height=5)

}
write.csv(moa.med,"synergy_forCorr.csv", row.names=FALSE)
