# analyzing IncuCyte results
library(plyr); library(dplyr); library(tidyr);library(ggplot2);library(colorspace)
setwd("~/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/MPNST-PDX-MT/drugViability")
dataPath <- "~/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/MPNST-PDX-MT/drugViability"

rel.conf <- read.csv("mpnst_combo_drug_response.csv")
musyc.scores <- read.csv("musyc_20250605.csv")

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
full.names <- c("MN-4","JH-2-079c", "JH-2-002")
short.names <- c("N-4","79c", "002")
for (i in 1:length(short.names)) {
  mean.conf[grepl(short.names[[i]],mean.conf$sample),]$sample <- full.names[[i]]
  musyc.scores[grepl(short.names[[i]],musyc.scores$sample),]$sample <- full.names[[i]]
}
prc2.loss <- c("JH-2-079c", "JH-2-002")
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
write.csv(musyc.scores, "Deconvolved_musyc_20250605.csv", row.names=FALSE)
musyc.scores <- read.csv("Deconvolved_musyc_20250605.csv")

# make scatter plot of a12 vs a21
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

musyc.scores$shortCombo <- paste0(tolower(substr(musyc.scores$drug1, 1,2)), "+",
                                  tolower(substr(musyc.scores$drug2, 1,2)))
ggplot(musyc.scores[musyc.scores$log_alpha12>0 & musyc.scores$log_alpha21>0,], 
       aes(x=log_alpha21, y=log_alpha12, shape=as.factor(time), color=drugCombo, size=R2)) +
  geom_point()+theme_classic() + facet_wrap(.~sample)+ scale_size_continuous(breaks=c(0.3,0.5,0.7)) +
  ggrepel::geom_text_repel(aes(label=shortCombo), box.padding = 0.5) + 
  labs(x="Drug 2's Effect on Potency of Drug 1", y="Drug 2's Effect on Potency of Drug 1", shape="Time (d)")
ggsave("musyc_potentCombos_v2.pdf",width=12, height=6)
#ggsave("musyc_potentCombos_classicTheme.pdf",width=12, height=6)

musyc.scores$meanLogAlpha <- NA
for (i in 1:nrow(musyc.scores)) {
  musyc.scores$meanLogAlpha[i] <- mean(c(musyc.scores$log_alpha12[i], musyc.scores$log_alpha21[i])) 
}
maxLogAlpha <- max(musyc.scores$meanLogAlpha)
ggplot(musyc.scores, aes(x=sample, y=reorder(drugCombo, meanLogAlpha), fill=meanLogAlpha)) + geom_tile(stat="identity") + 
  facet_wrap(.~paste0(time,"d"))+theme_classic() + scale_fill_gradient(low="grey",high="red",limits=c(0,8)) + 
  theme(axis.title = element_blank(), axis.text.x = element_text(angle=45, vjust=1, hjust=1)) + labs(fill="Mean\nLog|Alpha|")
ggsave("musyc_meanLogAlpha_heatmap_positive.pdf", width=4,height=4)

ggplot(musyc.scores, aes(x=sample, y=reorder(drugCombo, meanLogAlpha), fill=meanLogAlpha)) + geom_tile(stat="identity") + 
  facet_wrap(.~paste0(time,"d"))+theme_classic() + scale_fill_gradient2(low="blue",mid="grey",high="red",limits=c(-8,8)) + 
  theme(axis.title = element_blank(), axis.text.x = element_text(angle=45, vjust=1, hjust=1)) + labs(fill="Mean\nLog|Alpha|")
ggsave("musyc_meanLogAlpha_heatmap.pdf", width=4,height=4)

# try looking at Bliss scores
bliss <- data.frame()

dir.create("bliss")
blissFiles <- synapser::syncFromSynapse(entity='syn68639935', path='bliss')
blissSyn <- c("syn68685534","syn68685535", "syn68685536", "syn68685537","syn68685538","syn68685541")
for (i in blissSyn){
  test.bliss <- read.csv(synapser::synGet(i)$path)
  bliss <- rbind(bliss, test.bliss)
}
bliss$drugCombo <- paste0(bliss$drug1, "+", bliss$drug2)
bliss$sample <- bliss$PDX
bliss$time <- bliss$timePoint..hr./24
bliss.musyc <- merge(bliss, musyc.scores, by=c("drugCombo","sample","time"))
max.bliss <- plyr::ddply(bliss, .(drugCombo, sample, time), summarize,
                         maxBliss = max(Bliss_synergy))

ggplot(bliss.musyc, aes(x=meanLogAlpha, y=Bliss_synergy, shape=sample, color=drugCombo)) + #, size=R2)) +
  geom_point()+theme_minimal() + 
  #ggrepel::geom_text_repel(aes(label=paste0(time,"d"))) + 
  labs(x="MuSyC Potency", y="Bliss Synergy")
ggsave("musyc_meanLogAlpha_vs_bliss_synergy.pdf", width=6, height=6)

ggplot(bliss.musyc, aes(x=meanLogAlpha, y=Bliss_synergy, shape=sample, color=drugCombo)) + #, size=R2)) +
  geom_point()+theme_minimal() + facet_grid(paste0(time,"d")~sample)+
  #ggrepel::geom_text_repel(aes(label=paste0(time,"d"))) + 
  labs(x="MuSyC Potency", y="Bliss Synergy")
ggsave("musyc_meanLogAlpha_vs_bliss_synergy_Faceted.pdf", width=6, height=6)

cor.test(bliss.musyc$meanLogAlpha, bliss.musyc$Bliss_synergy)
cor.test(bliss.musyc$meanLogAlpha, bliss.musyc$Bliss_synergy, method="spearman")

bliss.musyc <- merge(max.bliss, musyc.scores, by=c("drugCombo","sample","time"))
ggplot(bliss.musyc, aes(x=meanLogAlpha, y=maxBliss, shape=sample, color=drugCombo)) + #, size=R2)) +
  geom_point()+theme_minimal() + 
  #ggrepel::geom_text_repel(aes(label=paste0(time,"d"))) + 
  labs(x="Mean MuSyC Potency", y="Max Bliss Synergy")
ggsave("musyc_meanLogAlpha_vs_max_bliss_synergy.pdf", width=6, height=6)

ggplot(bliss.musyc, aes(x=meanLogAlpha, y=maxBliss, shape=sample, color=drugCombo)) + #, size=R2)) +
  geom_point()+theme_minimal() + facet_grid(paste0(time,"d")~sample)+
  #ggrepel::geom_text_repel(aes(label=paste0(time,"d"))) + 
  labs(x="Mean MuSyC Potency", y="Max Bliss Synergy")
ggsave("musyc_meanLogAlpha_vs_max_bliss_synergy_Faceted.pdf", width=6, height=6)

cor.test(bliss.musyc$meanLogAlpha, bliss.musyc$maxBliss)
cor.test(bliss.musyc$meanLogAlpha, bliss.musyc$maxBliss, method="spearman")

# musyc.scores$universalCombo <- ""
# sampleCombos <- plyr::ddply(musyc.scores, .(drugCombo), summarize,
#                             samplesPos = sample[log_alpha12>0 & log_alpha21>0],
#                             )

ggplot(musyc.scores[musyc.scores$log_alpha12>0 & musyc.scores$log_alpha21>0,], 
       aes(x=log_alpha21, y=log_alpha12, shape=as.factor(time), color=drugCombo, size=R2)) +
  geom_point()+theme_classic() + facet_wrap(.~sample)+ scale_size_continuous(breaks=c(0.3,0.5,0.7)) +
  ggrepel::geom_text_repel(aes(label=shortCombo), box.padding = 0.5) + 
  labs(x="Drug 2's Effect on Potency of Drug 1", y="Drug 2's Effect on Potency of Drug 1", shape="Time (d)")
ggsave("musyc_potentCombos_v2.pdf",width=12, height=6)

ggplot(musyc.scores[musyc.scores$log_alpha12>0 & musyc.scores$log_alpha21>0 & musyc.scores$R2>=0.7,], 
       aes(x=log_alpha21, y=log_alpha12, shape=as.factor(time), color=drugCombo, size=R2)) +
  geom_point()+theme_minimal() + facet_wrap(.~sample)+
  ggrepel::geom_text_repel(aes(label=shortCombo), box.padding = 0.5) +
  labs(x="Drug 2's Effect on Potency of Drug 1", y="Drug 2's Effect on Potency of Drug 1", shape="Time (d)") #+
  #theme(legend.position="bottom"#, legend.direction="horizontal"
        #)
#ggsave("musyc_potentCombos_R20.7min_v3.pdf",width=10, height=4)
ggsave("musyc_potentCombos_R20.7min_v3.pdf",width=12, height=6)

# find max log alpha
musyc.scores[,c("mean_log_alpha","min_log_alpha","max_log_alpha")] <- NA
for (i in 1:nrow(musyc.scores)) {
  ci12 <- strsplit(musyc.scores$log_alpha12_ci[1],", ")[[1]]
  ci12 <- gsub("[[]","",ci12)
  ci12 <- gsub("[]]","",ci12)
  ci12 <- as.numeric(ci12)
  
  ci21 <- strsplit(musyc.scores$log_alpha21_ci[1],", ")[[1]]
  ci21 <- gsub("[[]","",ci21)
  ci21 <- gsub("[]]","",ci21)
  ci21 <- as.numeric(ci21)
  
  musyc.scores$mean_log_alpha[i] <- mean(c(musyc.scores$log_alpha12[i],
                                         musyc.scores$log_alpha21[i]))
  musyc.scores$max_log_alpha[i] <- max(c(ci12,ci21))
  musyc.scores$min_log_alpha[i] <- min(c(ci12,ci21))
}


results <- merge(mean.conf, musyc.scores, by=c("sample","time","drug1","drug2"), all.x=TRUE)

dir.create(paste0("curves_",Sys.Date()))
setwd(paste0("curves_",Sys.Date()))
library(drc) # curve source: answer by greenjune: https://stackoverflow.com/questions/36780357/plotting-dose-response-curves-with-ggplot2-and-drc
# Sara shared these 2 links: https://stackoverflow.com/questions/68209998/plot-drc-using-ggplot; https://forum.posit.co/t/extract-out-points-in-dose-response-curve-produced-by-drc/159433
results$timeD <- paste0(results$time,"d")
library(RColorBrewer)
cols1=c("#000000",brewer.pal(8,'Dark2'),brewer.pal(15-8,'Set2'),"mediumorchid1", "cornflowerblue")
drug.info <- list("palbociclib" = "CDK inhibitor", "ribociclib" = "CDK inhibitor",
                  "trabectedin" = "Chemotherapy", "Ifosfamide" = "DNA alkylating agent",
                  "decitabine" = "DNMT inhibitor", "vorinostat" = "HDAC inhibitor",
                  "mirdametinib" = "MEK inhibitor", "selumetinib" = "MEK inhibitor",
                  "trametinib" = "MEK inhibitor", "capmatinib" = "MET inhibitor",
                  "olaparib" = "PARP inhibitor", "RMC4630" = "SHP2 inhibitor",
                  "TNO155" = "SHP2 inhibitor", "doxorubicin" = "TOP inhibitor",
                  "irinotecan" = "TOP inhibitor", "verteporfin" = "YAP inhibitor",
                  "pexidartinib" = "KIT inhibitor")
names(cols1) <- c("DMSO", names(drug.info))
scales::show_col(cols1[drugs])
scales::show_col(cols1)
drugs[!(drugs %in% names(cols1))] # none missing
cols <- colorspace::darken(cols1, 0.3)
scales::show_col(cols)

##### just plot beta vs. logAlpha with R2
musyc.scores$Drugs <- paste0(musyc.scores$drug1,"+",musyc.scores$drug2)
musyc.scores$timeD <- paste0(musyc.scores$time,"d")
ggplot(musyc.scores, aes(x=beta, y=log_alpha12, color=R2, label=Drugs)) + 
  geom_point() + ggrepel::geom_text_repel() +
  facet_grid(timeD ~ sample) + theme_minimal()

ggplot(musyc.scores[musyc.scores$beta < 0 & musyc.scores$R2 > 0.8,], 
       aes(x=beta, y=log_alpha12, color=R2, label=Drugs)) + 
  geom_point() + ggrepel::geom_text_repel() +
  facet_grid(timeD ~ sample) + theme_minimal()

ggplot(musyc.scores[musyc.scores$beta < 0 & musyc.scores$R2 > 0.8,], 
       aes(x=beta, y=R2, label=Drugs)) + 
  geom_point() + ggrepel::geom_text_repel() +
  facet_grid(timeD ~ sample) + theme_minimal()

musyc.scores <- musyc.scores %>% tidyr::separate_wider_delim(beta_ci, ", ", names=c("beta_min", "beta_max"))
musyc.scores$beta_min <- substr(musyc.scores$beta_min,2,nchar(musyc.scores$beta_min))
musyc.scores$beta_max <- substr(musyc.scores$beta_max,1,nchar(musyc.scores$beta_max)-1)
class(musyc.scores$beta_min)
musyc.scores$beta_min <- as.numeric(musyc.scores$beta_min)
musyc.scores$beta_max <- as.numeric(musyc.scores$beta_max)
ggplot(musyc.scores[musyc.scores$beta < 0 & musyc.scores$R2 > 0.8,], 
       aes(x=Drugs, y=beta, fill=R2)) + 
  geom_bar(stat="identity") + 
  geom_errorbar(aes(ymin=beta_min, ymax=beta_max), width=0.2) +
  facet_grid(timeD ~ sample) + theme_minimal() + theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1))

ggplot(musyc.scores[musyc.scores$beta < 0 & musyc.scores$beta_max < 0 & musyc.scores$R2 > 0.8,], 
       aes(x=Drugs, y=beta, fill=R2)) + 
  geom_bar(stat="identity") + #scale_y_continuous(limits=c(-1,1))+
  geom_errorbar(aes(ymin=beta_min, ymax=beta_max)) +
  facet_grid(timeD ~ sample) + theme_minimal() + theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1))

musyc.scores$drugsShort <- paste0(substr(musyc.scores$drug1, 1, 3), "+", substr(musyc.scores$drug2, 1, 3))
plot.df <- musyc.scores[musyc.scores$beta < 0 & musyc.scores$beta_max < 0 & musyc.scores$beta_min<0 &
                          musyc.scores$beta_min < musyc.scores$beta & musyc.scores$beta_max > musyc.scores$beta &
                          musyc.scores$R2 > 0.7,]
drugOrder <- unique(plot.df[order(plot.df$beta),]$drugsShort)
ggplot(plot.df, aes(x=drugsShort, y=beta, fill=R2)) + 
  geom_bar(stat="identity") + #scale_y_continuous(limits=c(-2,0))+
  geom_errorbar(aes(ymin=beta_max, ymax=beta_min)) + 
  labs(y="% Change in Viability with Combination\nvs. Most Effective Single Agent") +
  facet_grid(timeD ~ sample) + theme_classic() + 
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1), axis.title.x=element_blank())
ggsave("musyc_beta_R2_0.7min.pdf", width=9, height=3)

plot.df <- musyc.scores[musyc.scores$beta < 0 & musyc.scores$beta_max < 0 & musyc.scores$beta_min<0 &
                          musyc.scores$beta_min < musyc.scores$beta & musyc.scores$beta_max > musyc.scores$beta,]
drugOrder <- unique(plot.df[order(plot.df$beta),]$drugsShort)
ggplot(plot.df, aes(x=drugsShort, y=beta, fill=R2)) + 
  geom_bar(stat="identity") + scale_y_continuous(limits=c(-2.2,0))+ scale_fill_continuous(limits=c(0,1)) +
  geom_errorbar(aes(ymin=beta_max, ymax=beta_min)) + 
  labs(y="% Change in Viability with Combination\nvs. Most Effective Single Agent") +
  facet_grid(timeD ~ sample) + theme_classic() + 
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1), axis.title.x=element_blank())
ggsave("musyc_beta_R2_yLim.pdf", width=9, height=3)

ggplot(plot.df, aes(x=drugsShort, y=beta, fill=R2)) + 
  geom_bar(stat="identity") + #scale_y_continuous(limits=c(-2.2,0))+ 
  scale_fill_continuous(limits=c(0,1)) +
  geom_errorbar(aes(ymin=beta_max, ymax=beta_min)) + 
  labs(y="% Change in Viability with Combination\nvs. Most Effective Single Agent") +
  facet_grid(timeD ~ sample) + theme_classic() + 
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1), axis.title.x=element_blank())
ggsave("musyc_beta_R2.pdf", width=9, height=3)

plot.df <- musyc.scores[musyc.scores$beta < 0 & musyc.scores$beta_max < 0 & musyc.scores$beta_min<0 &
                            musyc.scores$beta_min < musyc.scores$beta & musyc.scores$beta_max > musyc.scores$beta &
                            musyc.scores$sample == "JH-2-002",]
drugOrder <- unique(plot.df[order(plot.df$beta),]$drugsShort)
ggplot(plot.df, 
       aes(x=drugsShort, y=beta, fill=R2)) + scale_x_discrete(limits=drugOrder) +
  geom_bar(stat="identity") + labs(y="% Change in Viability with Combination\nvs. Most Effective Single Agent") +
  geom_errorbar(aes(ymin=beta_max, ymax=beta_min)) + ggtitle("JH-2-002") + scale_fill_continuous(limits=c(0,1)) +
  facet_grid(. ~ timeD) + theme_classic() + theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1),
                                                  plot.title=element_text(hjust=0.5, face="bold"), 
                                                  axis.title.x=element_blank())
ggsave("musyc_beta_R2_JH-2-002.pdf", width=6, height=3)

plot.df <- musyc.scores[musyc.scores$beta < 0 & musyc.scores$beta_max < 0 & musyc.scores$beta_min<0 &
                          musyc.scores$beta_min < musyc.scores$beta & musyc.scores$beta_max > musyc.scores$beta &
                          musyc.scores$sample == "JH-2-002" & musyc.scores$time == 2,]
drugOrder <- unique(plot.df[order(plot.df$beta),]$Drugs)
ggplot(plot.df, aes(x=Drugs, y=beta, fill=R2)) + scale_x_discrete(limits=drugOrder) +
  geom_bar(stat="identity") + scale_fill_continuous(limits=c(0,1))+ 
  labs(y="% Change in Viability with Combination\nvs. Most Effective Single Agent") +
  geom_errorbar(aes(ymin=beta_max, ymax=beta_min)) + ggtitle("JH-2-002: 2d") +
  theme_classic() + theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1),
                          plot.title=element_text(hjust=0.5, face="bold"), 
                          axis.title.x=element_blank())
ggsave("musyc_beta_R2_JH-2-002_2d.pdf", width=4, height=3)
ggsave("musyc_beta_R2_JH-2-002_2d_taller.pdf", width=4, height=4)

drugOrder <- unique(plot.df[order(plot.df$beta),]$drugsShort)
ggplot(plot.df, aes(x=drugsShort, y=beta, fill=R2)) + scale_x_discrete(limits=drugOrder) +
  geom_bar(stat="identity") + scale_fill_continuous(limits=c(0,1))+ 
  labs(y="% Change in Viability with Combination\nvs. Most Effective Single Agent") +
  geom_errorbar(aes(ymin=beta_max, ymax=beta_min)) + ggtitle("JH-2-002: 2d") +
  theme_classic() + theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1),
                          plot.title=element_text(hjust=0.5, face="bold"), 
                          axis.title.x=element_blank())
ggsave("musyc_beta_R2_JH-2-002_2d_shortNames.pdf", width=4, height=3)

hist(musyc.scores$beta)
hist(musyc.scores[abs(musyc.scores$beta) < 50,]$beta)
hist(musyc.scores[abs(musyc.scores$beta) < 20,]$beta)

##### ratio centric #####
results$conc12.ratio <- signif(as.numeric(results$drug1.conc/results$drug2.conc),3)
results$conc21.ratio <- signif(as.numeric(results$drug2.conc/results$drug1.conc),3)
results$combo.conc <- results$drug1.conc + results$drug2.conc
results$drugCombo <- paste0(results$drug1,"+",results$drug2)
library(plyr);library(dplyr);library(tidyr)
results <- results %>% tidyr::separate_wider_delim(beta_ci, ", ", names=c("beta_min", "beta_max"))
results$beta_min <- substr(results$beta_min,2,nchar(results$beta_min))
results$beta_max <- substr(results$beta_max,1,nchar(results$beta_max)-1)
class(results$beta_min)
results$beta_min <- as.numeric(results$beta_min)
results$beta_max <- as.numeric(results$beta_max)
results$potentialSynergy <- FALSE
results[!is.na(results$mean_log_alpha) & results$log_alpha12>0 & 
          results$log_alpha21>0,]$potentialSynergy <- TRUE
results$textFace <- "plain"
results[results$potentialSynergy,]$textFace <- "bold"

rel.conf$drugCombo <- paste0(rel.conf$drug1,"+",rel.conf$drug2)

all.p.df <- data.frame()
combos <- unique(results$drugCombo) # 17
for (c in combos) {
  combo.res <- results[results$drugCombo == c,]
  d1 <- unique(combo.res$drug1)
  d2 <- unique(combo.res$drug2)
  minConc <- min(combo.res$combo.conc)
  if (nrow(combo.res) > 0) {
    # generate plots
    temp.cols <- colorRampPalette(c(cols[[d2]], cols[[d1]]))(length(unique(combo.res$conc12.ratio)))
    mid.color <- temp.cols[ceiling(length(temp.cols)/2)]
    
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
    } else {
      combo.resp <- data.frame()
    }
    
    if (nrow(combo.resp) > 0) {
      if (any(combo.resp$p > 0.05)) {combo.resp[combo.resp$p > 0.05,]$textFace <- "plain"}
      conf.plot <- ggplot(combo.resp, aes(x=combo.conc, y=100*meanGROWTH, color=as.factor(conc12.ratio))) +
        facet_grid(timeD ~ sample) +
        geom_smooth(linetype="dashed", se=FALSE, method=drc::drm, method.args=list(fct=L.4(), se=FALSE)) +
        #scale_color_manual(values=c(scales::hue_pal()(2))) + 
        #scale_shape_manual(values=1:length(all.mpnst), breaks=all.mpnst)+
        scale_x_continuous(transform="log10") + geom_point() + 
        geom_errorbar(aes(ymin=100*(meanGROWTH-sdGROWTH), ymax=100*(meanGROWTH+sdGROWTH))) +
        #ggrepel::geom_text_repel(combo.res, aes(x=0,y=0,label=title))+
        ggtitle(paste0(d1," + ",d2))+scale_color_manual(values=temp.cols) +
        theme_classic(base_size=12) + labs(x="Concentration (uM)", y = "% Viability",
                                           color = paste("Ratio of",d1,"\nto",d2)) +
        theme(plot.title=element_text(hjust=0.5,face="bold"), 
              axis.text.x=element_text(angle=45, hjust=1, vjust=1))
      # conf.plotB <- conf.plot + geom_text(data=combo.resp, 
      #                                     mapping=aes(x=Inf, y=Inf, fontface=textFace,
      #                                                 label=paste0("Log|",as.character("\u03b1|=("),
      #                                                              signif(log_alpha12,2),",",signif(log_alpha21,2),
      #                                                              ")\n[",signif(min_log_alpha,2),",",signif(max_log_alpha,2),
      #                                                              "],\nR\u00b2=",signif(R2,2),
      #                                                              "\np=",signif(p,2))),
      #                                     hjust=1, vjust=1, show.legend=FALSE, 
      #                                     colour=mid.color)
                                          
      conf.plotB <- conf.plot + geom_text(data=combo.resp, 
                                          mapping=aes(x=Inf, y=Inf, fontface=textFace,
                                                      label=paste0("Log|",as.character("\u03b1|=("),
                                                                   signif(log_alpha12,2),",",signif(log_alpha21,2),
                                                                   "),\nR\u00b2=",signif(R2,2),
                                                                   ",\np=",signif(p,2))),
                                          hjust=1, vjust=1, show.legend=FALSE, 
                                          colour=mid.color)
      ggsave(paste0(d2,"_",d1,"_viability_log10_ratios_allMusyc_p.svg"),conf.plotB,width=9,height=9, device="svg") # was height 4
      ggsave(paste0(d2,"_",d1,"_viability_log10_ratios_allMusyc_p.png"),conf.plotB,width=9,height=9, device="png") # was height 4
    }
  }
}
write.csv(all.p.df, "tTests_viability_combo_lessThan_single.csv", row.names=FALSE)

#### compare musyc beta to DMEA ####
##### drug corr #####
synapser::synLogin()
drug.corr <- read.csv(synapser::synGet("syn65672266")$path)
drug.corr[drug.corr$DrugTreatment == "Trabectidin",]$DrugTreatment <- "Trabectedin"

shared.drugs <- drugs[tolower(drugs) %in% tolower(c(drug.corr$Drug, drug.corr$DrugTreatment))] # 9: mirdametinib, olaparib, TNO155, capmatinib, doxorubicin, palbociclib, vorinostat, irinotecan, trabectedin
shared.combos <- tidyr::expand_grid(shared.drugs, shared.drugs)
shared.combos$drugCombo <- paste0(shared.combos$shared.drugs...1,"+",shared.combos$shared.drugs...2)
shared.combos$drugCombo2 <- paste0(shared.combos$shared.drugs...2,"+",shared.combos$shared.drugs...1)
shared.combos <- unique(c(shared.combos$drugCombo, shared.combos$drugCombo2))
shared.combos <- shared.combos[shared.combos %in% combos] # 15 out of 17 tested combos

shared.corr <- drug.corr[tolower(drug.corr$Drug) %in% tolower(drugs) & 
                           tolower(drug.corr$DrugTreatment) %in% tolower(drugs),]
shared.corr$drugCombo <- NA
for (i in 1:nrow(shared.corr)) {
  d1 <- shared.corr$DrugTreatment[i]
  d2 <- shared.corr$Drug[i]
  
  # make lowercase unless TNO155
  d1 <- ifelse(d1=="TNO155",d1,tolower(d1))
  d2 <- ifelse(d2=="TNO155",d2,tolower(d2))
  
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
musyc.drug.corr <- cor.test(drug.corr.res$mean_log_alpha, drug.corr.res$Pearson.est)
musyc.drug.corr$estimate # 0.2528057
musyc.drug.corr$p.value # 7.345162e-17
musyc.drug.corrsp <- cor.test(drug.corr.res$mean_log_alpha, drug.corr.res$Pearson.est,method="spearman")
musyc.drug.corrsp$estimate # 0.2519774 
musyc.drug.corrsp$p.value # 9.323716e-17

musyc.drug.corr2 <- cor.test(drug.corr.res$mean_log_alpha, drug.corr.res$Spearman.est)
musyc.drug.corr2$estimate # 0.3478872 
musyc.drug.corr2$p.value # 2.096736e-31
musyc.drug.corr2sp <- cor.test(drug.corr.res$mean_log_alpha, drug.corr.res$Spearman.est,method="spearman")
musyc.drug.corr2sp$estimate # 0.3316385 
musyc.drug.corr2sp$p.value # 1.581662e-28

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
combo.cols <- unique(combo.cols) # 11

stats_pearson <- substitute(
  r == est * "," ~ ~"p" ~ "=" ~ p,
  list(
    est = format(as.numeric(musyc.drug.corr$estimate), digits = 2),
    p = format(musyc.drug.corr$p.value, digits = 2)
  )
)
ggplot(drug.corr.res, aes(x=mean_log_alpha, y=Pearson.est)) + # could set shape=sample but only sample is JH-2-002 
  geom_point(aes(color=drugCombo), show.legend =TRUE) + theme_minimal() + scale_color_manual(values=combo.cols, breaks=names(combo.cols),labels=names(combo.cols)) + 
  geom_smooth(method="lm",linetype="dashed",color="black") + ggplot2::geom_text(
    x = -Inf, y = Inf, vjust = "inward", hjust = "inward",
    colour = "black", parse = TRUE, 
    label = as.character(as.expression(stats_pearson)), size = 8
  ) +
  labs(x="MuSyC Synergy Score",
       y="DMEA Sensitivity Prediction")
ggsave("drugResults_musyc_meanLogAlpha_vs_dmea_PearsonEst.pdf",width=5,height=5)

stats_spearman <- substitute(
  rho == est * "," ~ ~"p" ~ "=" ~ p,
  list(
    est = format(as.numeric(musyc.drug.corrsp$estimate), digits = 2),
    p = format(musyc.drug.corrsp$p.value, digits = 2)
  )
)
ggplot(drug.corr.res, aes(x=mean_log_alpha, y=Pearson.est)) + # could set shape=sample but only sample is JH-2-002 
  geom_point(aes(color=drugCombo), show.legend =TRUE) + theme_minimal() + scale_color_manual(values=combo.cols, breaks=names(combo.cols),labels=names(combo.cols)) + 
  geom_smooth(method="lm",linetype="dashed",color="black") + ggplot2::geom_text(
    x = -Inf, y = Inf, vjust = "inward", hjust = "inward",
    colour = "black", parse = TRUE, 
    label = as.character(as.expression(stats_spearman)), size = 8
  ) +
  labs(x="MuSyC Synergy Score",
       y="DMEA Sensitivity Prediction")
ggsave("drugResults_musyc_meanLogAlpha_vs_dmea_PearsonEst_spearman.pdf",width=5,height=5)

shared.combos2 <- unique(drug.corr.res$drugCombo) # 10
ggsave("drugResults_musyc_meanLogAlpha_vs_dmea_PearsonEst.pdf",width=5,height=5)

stats_pearson <- substitute(
  r == est * "," ~ ~"p" ~ "=" ~ p,
  list(
    est = format(as.numeric(musyc.drug.corr2$estimate), digits = 2),
    p = format(musyc.drug.corr2$p.value, digits = 2)
  )
)
ggplot(drug.corr.res, aes(x=mean_log_alpha, y=Spearman.est)) + # could set shape=sample but only sample is JH-2-002 
  geom_point(aes(color=drugCombo), show.legend =TRUE) + theme_minimal() + scale_color_manual(values=combo.cols, breaks=names(combo.cols),labels=names(combo.cols)) + 
  geom_smooth(method="lm",linetype="dashed",color="black") + ggplot2::geom_text(
    x = -Inf, y = Inf, vjust = "inward", hjust = "inward",
    colour = "black", parse = TRUE, 
    label = as.character(as.expression(stats_pearson)), size = 8
  ) +
  labs(x="MuSyC Synergy Score",
       y="DMEA Sensitivity Prediction")
ggsave("drugResults_musyc_meanLogAlpha_vs_dmea_SpearmanEst.pdf",width=5,height=5)

stats_spearman <- substitute(
  rho == est * "," ~ ~"p" ~ "=" ~ p,
  list(
    est = format(as.numeric(musyc.drug.corr2sp$estimate), digits = 2),
    p = format(musyc.drug.corr2sp$p.value, digits = 2)
  )
)
ggplot(drug.corr.res, aes(x=mean_log_alpha, y=Spearman.est)) + # could set shape=sample but only sample is JH-2-002 
  geom_point(aes(color=drugCombo), show.legend =TRUE) + theme_minimal() + scale_color_manual(values=combo.cols, breaks=names(combo.cols),labels=names(combo.cols)) + 
  geom_smooth(method="lm",linetype="dashed",color="black") + ggplot2::geom_text(
    x = -Inf, y = Inf, vjust = "inward", hjust = "inward",
    colour = "black", parse = TRUE, 
    label = as.character(as.expression(stats_spearman)), size = 8
  ) +
  labs(x="MuSyC Synergy Score",
       y="DMEA Sensitivity Prediction")
ggsave("drugResults_musyc_meanLogAlpha_vs_dmea_SpearmanEst_spearman.pdf",width=5,height=5)
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

results2 <- results[,c("sample","drug1","drug2","mean_log_alpha","R2","timeD")]
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
moa.combos <- unique(results3$moaCombo) # 3

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

ggplot(shared.moa.res, aes(x=mean_log_alpha, y=NES, #shape=sample, # only sample is JH-2-002
                                                        color=moaComboShort)) +
  geom_point() + theme_minimal() + #scale_color_manual(values=combo.cols, breaks=names(combo.cols), labels=names(combo.cols)) +
  labs(x="MuSyC: % Change in Viability with Combination vs. Most Effective Single Agent",
       y="DMEA: Predicted Drug Sensitivity")
musyc.drug.corr <- cor.test(shared.moa.res$mean_log_alpha, shared.moa.res$NES)
musyc.drug.corr$estimate # -0.07149301 
musyc.drug.corr$p.value # 0.08647383
musyc.drug.corrsp <- cor.test(shared.moa.res$mean_log_alpha, shared.moa.res$NES ,method="spearman")
musyc.drug.corrsp$estimate # -0.09187795 
musyc.drug.corrsp$p.value # 0.02745786

ggplot(shared.moa.res[shared.moa.res$sig,], aes(x=mean_log_alpha, y=NES, #shape=sample, # only sample is JH-2-002
                           color=moaComboShort)) +
  geom_point() + theme_minimal() + #scale_color_manual(values=combo.cols, breaks=names(combo.cols), labels=names(combo.cols)) +
  labs(x="MuSyC: % Change in Viability with Combination vs. Most Effective Single Agent",
       y="DMEA: Predicted Drug Sensitivity")
musyc.drug.corr <- cor.test(shared.moa.res[shared.moa.res$sig,]$mean_log_alpha, shared.moa.res[shared.moa.res$sig,]$NES)
musyc.drug.corr$estimate # 0.09917356
musyc.drug.corr$p.value # 0.09298546
musyc.drug.corrsp <- cor.test(shared.moa.res[shared.moa.res$sig,]$mean_log_alpha, shared.moa.res[shared.moa.res$sig,]$NES ,method="spearman")
musyc.drug.corrsp$estimate # 0
musyc.drug.corrsp$p.value # 1

ggplot(shared.moa.res2, aes(y=mean_log_alpha, x=NES)) + 
  geom_point(aes(shape=moaComboShort, color=sig), size=3) + 
  theme_minimal() + ggrepel::geom_text_repel(data=subset(shared.moa.res2,sig),
                                             aes(label=paste(as.character(Time), "RNA,\n",as.character(timeD),"viability"))) + 
  labs(y="MuSyC Synergy Score", shape = "Drug Combination", 
       x="DMEA Sensitivity Score", color = "Significant\nDMEA Result") + 
  scale_color_manual(values=c("red","grey"), breaks=c(TRUE, FALSE))
ggsave("MOA_results_musyc_vs_dmea.pdf", width=7, height=4)

ggplot(shared.moa.res2, aes(y=mean_log_alpha, x=NES)) + 
  geom_point(aes(shape=moaComboShort, color=sig), size=3) + 
  theme_minimal() + 
  labs(y="MuSyC Synergy Score", shape = "Drug Combination", 
       x="DMEA Sensitivity Score", color = "Significant\nDMEA Result") + 
  scale_color_manual(values=c("red","grey"), breaks=c(TRUE, FALSE))
ggsave("MOA_results_musyc_vs_dmea_noLabels.pdf", width=7, height=4)

ggplot(shared.moa.res2, aes(y=mean_log_alpha, x=NES)) + 
  geom_point(aes(shape=moaComboShort, color=sig), size=3) + 
  geom_smooth(method="lm",se=FALSE,formula = y ~ x + I(x^2), color="darkgrey", linetype="dashed") +
  theme_minimal() + 
  labs(y="MuSyC Synergy Score", shape = "Drug Combination", 
       x="DMEA Sensitivity Score", color = "Significant\nDMEA Result") +
  scale_color_manual(values=c("red","grey"), breaks=c(TRUE, FALSE))
ggsave("MOA_results_musyc_vs_dmea_noLabels_wQuadraticFit.pdf", width=7, height=4)

ggplot(shared.moa.res2, aes(y=mean_log_alpha, x=NES)) + 
  geom_point(aes(shape=moaComboShort, color=sig), size=3) + 
  geom_smooth(method="lm",se=FALSE,formula = y ~ x + I(x^2), color="darkgrey", linetype="dashed") +
  theme_minimal() + ggrepel::geom_text_repel(data=subset(shared.moa.res2,sig),
                                             aes(label=paste(as.character(Time), "RNA,\n",as.character(timeD),"viability"))) + 
  labs(y="MuSyC Synergy Score", shape = "Drug Combination", 
       x="DMEA Sensitivity Score", color = "Significant\nDMEA Result") +
  scale_color_manual(values=c("red","grey"), breaks=c(TRUE, FALSE))
ggsave("MOA_results_musyc_vs_dmea_wQuadraticFit.pdf", width=7, height=4)
write.csv(shared.moa.res2,"MOA_results_musyc_vs_dmea.csv", row.names=FALSE)

ggplot(shared.moa.res2[shared.moa.res2$mean_log_alpha>0,], aes(y=mean_log_alpha, x=NES)) + 
  geom_point(aes(shape=moaComboShort, color=sig), size=3) + 
  theme_minimal() + ggrepel::geom_text_repel(data=subset(shared.moa.res2[shared.moa.res2$mean_log_alpha>0,],sig),
                                             aes(label=paste(as.character(Time), "RNA,\n",as.character(timeD),"viability"))) + 
  labs(y="MuSyC Synergy Score", shape = "Drug Combination", 
       x="DMEA Sensitivity Score", color = "Significant\nDMEA Result") + 
  scale_color_manual(values=c("red","grey"), breaks=c(TRUE, FALSE))
ggsave("MOA_results_musyc_vs_dmea_positiveAlpha.pdf", width=7, height=4)

ggplot(shared.moa.res2[shared.moa.res2$mean_log_alpha>0,], aes(y=mean_log_alpha, x=NES)) + 
  geom_point(aes(shape=moaComboShort, color=sig), size=3) + 
  geom_smooth(method="lm",se=FALSE,formula = y ~ x + I(x^2), color="darkgrey", linetype="dashed") +
  theme_minimal() + ggrepel::geom_text_repel(data=subset(shared.moa.res2[shared.moa.res2$mean_log_alpha>0,],sig),
                                             aes(label=paste(as.character(Time), "RNA,\n",as.character(timeD),"viability"))) + 
  labs(y="MuSyC Synergy Score", shape = "Drug Combination", 
       x="DMEA Sensitivity Score", color = "Significant\nDMEA Result") +
  scale_color_manual(values=c("red","grey"), breaks=c(TRUE, FALSE))
ggsave("MOA_results_musyc_vs_dmea_wQuadraticFit_positiveAlpha.pdf", width=7, height=4)

ggplot(shared.moa.res2[shared.moa.res2$mean_log_alpha>0,], aes(y=mean_log_alpha, x=NES)) + 
  geom_point(aes(shape=moaComboShort, color=sig), size=3) + 
  geom_smooth(method="lm",se=FALSE,formula = y ~ x + I(x^2), color="darkgrey", linetype="dashed") +
  theme_minimal() + 
  labs(y="MuSyC Synergy Score", shape = "Drug Combination", 
       x="DMEA Sensitivity Score", color = "Significant\nDMEA Result") +
  scale_color_manual(values=c("red","grey"), breaks=c(TRUE, FALSE))
ggsave("MOA_results_musyc_vs_dmea_wQuadraticFit_positiveAlpha_noLabels.pdf", width=7, height=4)

# are mean_log_alpha values higher in significant DMEA results?
b.test <- t.test(shared.moa.res2[shared.moa.res2$sig,]$mean_log_alpha, shared.moa.res2[!shared.moa.res2$sig,]$mean_log_alpha, alternative="greater")
b.test$p.value # 0.03112967
b.test$estimate # mean of x mean of y 
# 2.8980339 0.4124156 
length(shared.moa.res2[shared.moa.res2$sig,]$mean_log_alpha) # 12
length(shared.moa.res2[!shared.moa.res2$sig,]$mean_log_alpha) # 12

b.test <- t.test(shared.moa.res2[shared.moa.res2$sig & shared.moa.res2$mean_log_alpha>0,]$mean_log_alpha, 
                 shared.moa.res2[!shared.moa.res2$sig & shared.moa.res2$mean_log_alpha>0,]$mean_log_alpha, alternative="greater")
b.test$p.value # 0.06034868
b.test$estimate # mean of x mean of y 
# 4.179305  3.074639 
length(shared.moa.res2[shared.moa.res2$sig & abs(shared.moa.res2$mean_log_alpha)<2,]$mean_log_alpha) # 3
length(shared.moa.res2[!shared.moa.res2$sig & abs(shared.moa.res2$mean_log_alpha)<2,]$mean_log_alpha) # 1