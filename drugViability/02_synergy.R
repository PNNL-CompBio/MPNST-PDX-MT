# compile synergy results
# author: Belinda B. Garana
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

cor.test(bliss.musyc$meanLogAlpha, bliss.musyc$maxBliss) # r = 0.011, p=0.89
cor.test(bliss.musyc$meanLogAlpha, bliss.musyc$maxBliss, method="spearman") # rho = 0.039, p = 0.6151

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

cor.test(bliss.musyc.tested$meanLogAlpha, bliss.musyc.tested$maxBliss) # r=0.066, p=0.3994
cor.test(bliss.musyc.tested$meanLogAlpha, bliss.musyc.tested$maxBliss, method="spearman") # rho = 0.1481161, p=0.0561

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