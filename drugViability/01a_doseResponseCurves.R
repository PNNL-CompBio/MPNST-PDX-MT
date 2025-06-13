# analyzing IncuCyte results
library(plyr); library(dplyr); library(tidyr);library(ggplot2)
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
results <- merge(mean.conf, musyc.scores, by=c("sample","time","drug1","drug2"), all.x=TRUE)
results$title <- NA
results[!is.na(results$beta) & results$beta > 0,]$title <- 
  paste0(signif(results[!is.na(results$beta) & results$beta > 0,]$beta,2), 
         "% increase (", expression("R"^2),"=",
         signif(results[!is.na(results$beta) & results$beta > 0,]$R2,2),")")
# results[!is.na(results$beta) & results$beta > 0,]$title <- paste0(signif(results[!is.na(results$beta) & results$beta > 0,]$beta,2), 
#                                            "% increase in efficacy with combination (R2=",
#                                            signif(results[!is.na(results$beta) & results$beta > 0,]$R2,2),")")
results[!is.na(results$beta) & results$beta < 0,]$title <- 
  paste0(signif(results[!is.na(results$beta) & results$beta < 0,]$beta,2), 
         "% decrease (", expression("R"^2),"=",
         signif(results[!is.na(results$beta) & results$beta < 0,]$R2,2),")") # in viability with combo vs. most efficacious single agent


dir.create("curves")
setwd("curves")
library(drc) # curve source: answer by greenjune: https://stackoverflow.com/questions/36780357/plotting-dose-response-curves-with-ggplot2-and-drc
# Sara shared these 2 links: https://stackoverflow.com/questions/68209998/plot-drc-using-ggplot; https://forum.posit.co/t/extract-out-points-in-dose-response-curve-produced-by-drc/159433
results$timeD <- paste0(results$time,"d")
drugs1 <- unique(results$drug1)
#### individual drugs or combos only ####
##### drug1 centric #####
for (d1 in drugs1) {
  d1.res <- results[results$drug1 == d1,]
  drugs2 <- unique(d1.res$drug2)
  doses1 <- unique(d1.res$drug1.conc)

  for (dose in doses1) {
    d1.dose.res <- d1.res[d1.res$drug1.conc == dose,]
    
    for (d2 in drugs2) {
      combo.res <- d1.dose.res[d1.dose.res$drug2 == d2,]

      if (nrow(combo.res) > 0 & length(unique(combo.res$drug2.conc))>2) {
        # if (any(combo.res$drug2.conc == 0)) {
        #   combo.res[combo.res$drug2.conc == 0,]$drug2.conc <- 1E-300
        # }
        # generate plots
        conf.plot <- ggplot(combo.res,
                            aes(x=drug2.conc, y=100*meanGROWTH, color=PRC2, shape=sample)) +
          facet_wrap(.~timeD) +
          geom_smooth(linetype="dashed", se=FALSE, method=drc::drm, method.args=list(fct=L.4(), se=FALSE)) +
          #scale_color_manual(values=c(scales::hue_pal()(2))) + 
          #scale_shape_manual(values=1:length(all.mpnst), breaks=all.mpnst)+
          scale_x_continuous(transform="log10") + geom_point() +
          geom_errorbar(aes(ymin=meanGROWTH-sdGROWTH, ymax=meanGROWTH+sdGROWTH), width=0.2) +
          #ggrepel::geom_text_repel(combo.res, aes(x=0,y=0,label=title))+
          ggtitle(paste0(d2," + ",signif(dose,2),"uM ",d1))+
          theme_classic(base_size=12) + labs(x="Concentration (uM)", y = "% Viability",
                                             color = "PRC2 Status", shape = "MPNST") +
          theme(plot.title=element_text(hjust=0.5,face="bold"))
        ggsave(paste0(d1,"-",signif(dose,2),"uM_",d2,"_viability_log10.pdf"),conf.plot,width=6,height=4)
      }
    }
  }
}

##### drug2 centric #####
drugs2 <- unique(results$drug2)
for (d2 in drugs2) {
  d2.res <- results[results$drug2 == d2,]
  drugs1 <- unique(d1.res$drug1)
  doses2 <- unique(d1.res$drug2.conc)
  
  for (dose in doses2) {
    d2.dose.res <- d2.res[d2.res$drug2.conc == dose,]
    
    for (d1 in drugs1) {
      combo.res <- d2.dose.res[d2.dose.res$drug1 == d1,]
      
      if (nrow(combo.res) > 0 & length(unique(combo.res$drug1.conc))>2) {
        # if (any(combo.res$drug1.conc == 0)) {
        #   combo.res[combo.res$drug1.conc == 0,]$drug1.conc <- 1E-300
        # }
        # generate plots
        conf.plot <- ggplot(combo.res,
                            aes(x=drug1.conc, y=100*meanGROWTH, color=PRC2, shape=sample)) +
          facet_wrap(.~timeD) +
          geom_smooth(linetype="dashed", se=FALSE, method=drc::drm, method.args=list(fct=L.4(), se=FALSE)) +
          #scale_color_manual(values=c(scales::hue_pal()(2))) + 
          #scale_shape_manual(values=1:length(all.mpnst), breaks=all.mpnst)+
          scale_x_continuous(transform="log10") + geom_point() +
          geom_errorbar(aes(ymin=meanGROWTH-sdGROWTH, ymax=meanGROWTH+sdGROWTH), width=0.2) +
          #ggrepel::geom_text_repel(combo.res, aes(x=0,y=0,label=title))+
          ggtitle(paste0(d2," + ",signif(dose,2),"uM ",d1))+
          theme_classic(base_size=12) + labs(x="Concentration (uM)", y = "% Viability",
                                             color = "PRC2 Status", shape = "MPNST") +
          theme(plot.title=element_text(hjust=0.5,face="bold"))
        ggsave(paste0(d2,"-",signif(dose,2),"uM_",d1,"_viability_log10.pdf"),conf.plot,width=6,height=4) 
      }
    }
  }
}

##### ratio centric #####
results$conc12.ratio <- as.numeric(results$drug1.conc/results$drug2.conc)
for (d1 in drugs1) {
  d1.res <- results[results$drug1 == d1,]
  drugs2 <- unique(d1.res$drug2)
  doses1 <- unique(d1.res$conc12.ratio)
  doses1 <- doses1[doses1 != "Inf"]
  
  for (dose in doses1) {
    d1.dose.res <- d1.res[d1.res$conc12.ratio == dose,]
    
    for (d2 in drugs2) {
      combo.res <- d1.dose.res[d1.dose.res$drug2 == d2,]
      
      if (nrow(combo.res) > 0) {
        # if (any(combo.res$drug1.conc == 0)) {
        #   combo.res[combo.res$drug1.conc == 0,]$drug1.conc <- 1E-300
        # }
        # if (any(combo.res$drug2.conc == 0)) {
        #   combo.res[combo.res$drug2.conc == 0,]$drug2.conc <- 1E-300
        # }
        # generate plots
        if (length(unique(combo.res$drug2.conc)) > 2) {
          conf.plot <- ggplot(combo.res,
                              aes(x=drug2.conc, y=100*meanGROWTH, color=PRC2, shape=sample)) +
            facet_wrap(.~timeD) +
            geom_smooth(linetype="dashed", se=FALSE, method=drc::drm, method.args=list(fct=L.4(), se=FALSE)) +
            #scale_color_manual(values=c(scales::hue_pal()(2))) + 
            #scale_shape_manual(values=1:length(all.mpnst), breaks=all.mpnst)+
            scale_x_continuous(transform="log10") + geom_point() +
            geom_errorbar(aes(ymin=meanGROWTH-sdGROWTH, ymax=meanGROWTH+sdGROWTH), width=0.2) +
            #ggrepel::geom_text_repel(combo.res, aes(x=0,y=0,label=title))+
            ggtitle(paste0(d1," + ",d2," (",signif(dose,2)," ratio)"))+
            theme_classic(base_size=12) + labs(x=paste(d2, "Concentration (uM)"), y = "% Viability",
                                               color = "PRC2 Status", shape = "MPNST") +
            theme(plot.title=element_text(hjust=0.5,face="bold"))
          ggsave(paste0(d1,"_",d2,"_",signif(dose,2),"concRatio_viability_log10Drug2.pdf"),conf.plot,width=6,height=4) 
        }
        
        if (length(unique(combo.res$drug1.conc)) > 2) {
          conf.plot <- ggplot(combo.res,
                              aes(x=drug1.conc, y=100*meanGROWTH, color=PRC2, shape=sample)) +
            facet_wrap(.~timeD) +
            geom_smooth(linetype="dashed", se=FALSE, method=drc::drm, method.args=list(fct=L.4(), se=FALSE)) +
            #scale_color_manual(values=c(scales::hue_pal()(2))) + 
            #scale_shape_manual(values=1:length(all.mpnst), breaks=all.mpnst)+
            scale_x_continuous(transform="log10") + geom_point() +
            geom_errorbar(aes(ymin=meanGROWTH-sdGROWTH, ymax=meanGROWTH+sdGROWTH), width=0.2) +
            #ggrepel::geom_text_repel(combo.res, aes(x=0,y=0,label=title))+
            ggtitle(paste0(d1," + ",d2," (",signif(dose,2)," ratio)"))+
            theme_classic(base_size=12) + labs(x=paste(d2, "Concentration (uM)"), y = "% Viability",
                                               color = "PRC2 Status", shape = "MPNST") +
            theme(plot.title=element_text(hjust=0.5,face="bold"))
          ggsave(paste0(d1,"_",d2,"_",signif(dose,2),"concRatio_viability_log10Drug1.pdf"),conf.plot,width=6,height=4) 
        }
      }
    }
  }
}

#### plot single and combo on same axes ####
dir.create("2025-06-11")
setwd("2025-06-11")
library(RColorBrewer)
cols=c("black",brewer.pal(8,'Dark2'),brewer.pal(15-8,'Set2'),"mediumorchid1", "cornflowerblue")
drug.info <- list("palbociclib" = "CDK inhibitor", "ribociclib" = "CDK inhibitor",
                  "trabectedin" = "Chemotherapy", "Ifosfamide" = "DNA alkylating agent",
                  "decitabine" = "DNMT inhibitor", "vorinostat" = "HDAC inhibitor",
                  "mirdametinib" = "MEK inhibitor", "selumetinib" = "MEK inhibitor",
                  "trametinib" = "MEK inhibitor", "capmatinib" = "MET inhibitor",
                  "olaparib" = "PARP inhibitor", "RMC4630" = "SHP2 inhibitor",
                  "TNO155" = "SHP2 inhibitor", "doxorubicin" = "TOP inhibitor",
                  "irinotecan" = "TOP inhibitor", "verteporfin" = "YAP inhibitor",
                  "pexidartinib" = "KIT inhibitor")
names(cols) <- c("DMSO", names(drug.info))
scales::show_col(cols[drugs])
scales::show_col(cols)
drugs[!(drugs %in% names(cols))] # none missing
results$drugs1 <- paste0(results$drug1," (",signif(results$drug1.conc,2),"uM)+",results$drug2)
results$drugs2 <- paste0(results$drug2," (",signif(results$drug2.conc,2),"uM)+",results$drug1)
results$drug <- paste0(results$drug1,"+",results$drug2," (", signif(results$conc12.ratio,2)," ratio)")
results[results$conc12.ratio == 0,]$drug <- results[results$conc12.ratio == 0,]$drug2
results[results$conc12.ratio == "Inf",]$drug <- results[results$conc12.ratio == "Inf",]$drug1
results[results$conc12.ratio == 0,]$drugs1 <- results[results$conc12.ratio == 0,]$drug2
results[results$conc12.ratio == "Inf",]$drugs1 <- results[results$conc12.ratio == "Inf",]$drug1
results[results$conc12.ratio == 0,]$drugs2 <- results[results$conc12.ratio == 0,]$drug2
results[results$conc12.ratio == "Inf",]$drugs2 <- results[results$conc12.ratio == "Inf",]$drug1

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
##### drug1 centric #####
#remotes::install_github("jmw86069/colorjam");
library(colorjam)
library(ggplot2)
library(drc)
drugs1 <- unique(results$drug1)
for (d1 in drugs1) {
  d1.res <- results[results$drug1 == d1,]
  drugs2 <- unique(d1.res$drug2)
  doses1 <- unique(d1.res$drug1.conc)
  
  for (d2 in drugs2) {
    combo.res <- d1.res[d1.res$drug2 == d2,]
    
    if (nrow(combo.res) > 0 & length(unique(combo.res$drug2.conc))>2) {
      # if (any(combo.res$drug2.conc == 0)) {
      #   combo.res[combo.res$drug2.conc == 0,]$drug2.conc <- 1E-300
      # }
      # generate plots
      temp.cols <- colorRampPalette(c(cols[[d2]], 
                                      colorjam::blend_colors(c(cols[[d1]],cols[[d2]]))))(
                                        length(unique(combo.res$drug1.conc)))
      conf.plot <- ggplot(combo.res,
                          aes(x=drug2.conc, y=100*meanGROWTH, color=as.factor(signif(drug1.conc,2)),
                          group=as.factor(signif(drug1.conc,2)))) +
        facet_grid(timeD ~ sample) +
        geom_smooth(linetype="dashed", se=FALSE, method=drc::drm, method.args=list(fct=L.4(), se=FALSE)) +
        #scale_color_manual(values=c(scales::hue_pal()(2))) + 
        #scale_shape_manual(values=1:length(all.mpnst), breaks=all.mpnst)+
        scale_x_continuous(transform="log10") + geom_point() + 
        geom_errorbar(aes(ymin=100*(meanGROWTH-sdGROWTH), ymax=100*(meanGROWTH+sdGROWTH))) +
        #ggrepel::geom_text_repel(combo.res, aes(x=0,y=0,label=title))+
        ggtitle(paste0(d2," + ",d1))+scale_color_manual(values=temp.cols) +
        theme_classic(base_size=12) + labs(x=paste(d2, "concentration (uM)"), y = "% Viability",
                                           color = paste0(d1, "\nconcentration (uM)")) +
        theme(plot.title=element_text(hjust=0.5,face="bold"))
      ggsave(paste0(d1,"_",d2,"_viability_log10.pdf"),conf.plot,width=8,height=4)

      temp.cols <- colorRampPalette(c(cols[[d2]],
                                      "gray88"))(
                                        length(unique(combo.res$drug1.conc)))
      conf.plot <- ggplot(combo.res,
                          aes(x=drug2.conc, y=100*meanGROWTH, color=as.factor(signif(drug1.conc,2)),
                          group=as.factor(signif(drug1.conc,2)))) +
        facet_grid(timeD ~ sample) +
        geom_smooth(linetype="dashed", se=FALSE, method=drc::drm, method.args=list(fct=L.4(), se=FALSE)) +
        #scale_color_manual(values=c(scales::hue_pal()(2))) +
        #scale_shape_manual(values=1:length(all.mpnst), breaks=all.mpnst)+
        scale_x_continuous(transform="log10") + geom_point() +
        geom_errorbar(aes(ymin=100*(meanGROWTH-sdGROWTH), ymax=100*(meanGROWTH+sdGROWTH))) +
        #ggrepel::geom_text_repel(combo.res, aes(x=0,y=0,label=title))+
        ggtitle(paste0(d2," + ",d1))+scale_color_manual(values=temp.cols) +
        theme_classic(base_size=12) + labs(x=paste(d2, "concentration (uM)"), y = "% Viability",
                                           color = paste0(d1, "\nconcentration (uM)")) +
        theme(plot.title=element_text(hjust=0.5,face="bold"))
      ggsave(paste0(d1,"_",d2,"_viability_log10_gray88_v2.pdf"),conf.plot,width=8,height=4)
      
      temp.cols2 <- colorRampPalette(c(cols[[d2]], 
                                      "gray88"))(3)
      drug1.conc.vals <- c(0,median(combo.res$drug1.conc), max(combo.res$drug1.conc))
      names(temp.cols2) <- signif(drug1.conc.vals,2)
      conf.plot <- ggplot(combo.res[combo.res$drug1.conc %in% drug1.conc.vals,],
                          aes(x=drug2.conc, y=100*meanGROWTH, color=as.factor(signif(drug1.conc,2)), 
                              group=as.factor(signif(drug1.conc,2)))) +
        facet_grid(timeD ~ sample) +
        geom_smooth(linetype="dashed", se=FALSE, method=drc::drm, method.args=list(fct=L.4(), se=FALSE)) +
        #scale_color_manual(values=c(scales::hue_pal()(2))) + 
        #scale_shape_manual(values=1:length(all.mpnst), breaks=all.mpnst)+
        scale_x_continuous(transform="log10") + geom_point() + 
        geom_errorbar(aes(ymin=100*(meanGROWTH-sdGROWTH), ymax=100*(meanGROWTH+sdGROWTH))) +
        #ggrepel::geom_text_repel(combo.res, aes(x=0,y=0,label=title))+
        ggtitle(paste0(d2," + ",d1))+scale_color_manual(values=temp.cols2, breaks=names(temp.cols2)) +
        theme_classic(base_size=12) + labs(x=paste(d2, "concentration (uM)"), y = "% Viability",
                                           color = paste0(d1, "\nconcentration (uM)")) +
        theme(plot.title=element_text(hjust=0.5,face="bold"))
      ggsave(paste0(d1,"_",d2,"_viability_log10_gray88_v3.pdf"),conf.plot,width=8,height=4)
      
      # conf.plot <- ggplot(combo.res,
      #                     aes(x=drug2.conc, y=100*meanGROWTH, color=as.factor(signif(drug1.conc,2)))) +
      #   facet_grid(timeD ~ sample) +
      #   geom_smooth(linetype="dashed", se=FALSE, method=drc::drm, method.args=list(fct=L.4(), se=FALSE)) +
      #   #scale_color_manual(values=c(scales::hue_pal()(2))) + 
      #   #scale_shape_manual(values=1:length(all.mpnst), breaks=all.mpnst)+
      #   scale_x_continuous(transform="log10") + geom_point() + 
      #   geom_errorbar(aes(ymin=meanGROWTH-sdGROWTH, ymax=meanGROWTH+sdGROWTH)) +
      #   #ggrepel::geom_text_repel(combo.res, aes(x=0,y=0,label=title))+
      #   ggtitle(paste0(d2," + ",d1))+scale_color_manual(values=temp.cols) +
      #   theme_classic(base_size=12) + labs(x=paste(d2, "concentration (uM)"), y = "% Viability",
      #                                      color = paste0(d1, "\nconcentration (uM)")) +
      #   theme(plot.title=element_text(hjust=0.5,face="bold"))
      # ggsave(paste0(d1,"_",d2,"_viability_log10_gray88_v2.pdf"),conf.plot,width=8,height=4)
    }
  }
}


##### drug2 centric #####
drugs2 <- unique(results$drug2)
for (d2 in drugs2) {
  d2.res <- results[results$drug2 == d2,]
  drugs1 <- unique(d2.res$drug1)
  doses2 <- unique(d2.res$drug2.conc)
  
  for (d1 in drugs1) {
    combo.res <- d2.res[d2.res$drug1 == d1,]
    
    if (nrow(combo.res) > 0 & length(unique(combo.res$drug1.conc))>2) {
      # generate plots
      temp.cols <- colorRampPalette(c(cols[[d1]], 
                                      colorjam::blend_colors(c(cols[[d2]],cols[[d1]]))))(
                                        length(unique(combo.res$drug2.conc)))
      conf.plot <- ggplot(combo.res,
                          aes(x=drug1.conc, y=100*meanGROWTH, color=as.factor(signif(drug2.conc,2)),
                              group=as.factor(signif(drug2.conc,2)))) +
        facet_grid(timeD ~ sample) +
        geom_smooth(linetype="dashed", se=FALSE, method=drc::drm, method.args=list(fct=L.4(), se=FALSE)) +
        #scale_color_manual(values=c(scales::hue_pal()(2))) + 
        #scale_shape_manual(values=1:length(all.mpnst), breaks=all.mpnst)+
        scale_x_continuous(transform="log10") + geom_point() + 
        geom_errorbar(aes(ymin=100*(meanGROWTH-sdGROWTH), ymax=100*(meanGROWTH+sdGROWTH))) +
        #ggrepel::geom_text_repel(combo.res, aes(x=0,y=0,label=title))+
        ggtitle(paste0(d1," + ",d2))+scale_color_manual(values=temp.cols) +
        theme_classic(base_size=12) + labs(x=paste(d1, "concentration (uM)"), y = "% Viability",
                                           color = paste0(d2, "\nconcentration (uM)")) +
        theme(plot.title=element_text(hjust=0.5,face="bold"))
      ggsave(paste0(d2,"_",d1,"_viability_log10.pdf"),conf.plot,width=8,height=4)
      
      temp.cols <- colorRampPalette(c(cols[[d1]],
                                      "gray88"))(
                                        length(unique(combo.res$drug2.conc)))
      conf.plot <- ggplot(combo.res,
                          aes(x=drug1.conc, y=100*meanGROWTH, color=as.factor(signif(drug2.conc,2)),
                              group=as.factor(signif(drug2.conc,2)))) +
        facet_grid(timeD ~ sample) +
        geom_smooth(linetype="dashed", se=FALSE, method=drc::drm, method.args=list(fct=L.4(), se=FALSE)) +
        #scale_color_manual(values=c(scales::hue_pal()(2))) +
        #scale_shape_manual(values=1:length(all.mpnst), breaks=all.mpnst)+
        scale_x_continuous(transform="log10") + geom_point() +
        geom_errorbar(aes(ymin=100*(meanGROWTH-sdGROWTH), ymax=100*(meanGROWTH+sdGROWTH))) +
        #ggrepel::geom_text_repel(combo.res, aes(x=0,y=0,label=title))+
        ggtitle(paste0(d1," + ",d2))+scale_color_manual(values=temp.cols) +
        theme_classic(base_size=12) + labs(x=paste(d1, "concentration (uM)"), y = "% Viability",
                                           color = paste0(d2, "\nconcentration (uM)")) +
        theme(plot.title=element_text(hjust=0.5,face="bold"))
      ggsave(paste0(d2,"_",d1,"_viability_log10_gray88.pdf"),conf.plot,width=8,height=4)
    }
  }
}

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
results[!is.na(results$beta) & results$beta<0 & results$beta_max<0 &
          results$beta<results$beta_max & 
          results$beta_min<results$beta,]$potentialSynergy <- TRUE
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
    conf.plot <- ggplot(combo.res, aes(x=combo.conc, y=100*meanGROWTH, color=as.factor(conc12.ratio))) +
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
    ggsave(paste0(d2,"_",d1,"_viability_log10_ratios.pdf"),conf.plot,width=9,height=4)
    
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
    
    if (any(combo.res$potentialSynergy)) {
      conf.plotB <- conf.plot + geom_text(data=subset(combo.res, potentialSynergy), 
                                          mapping=aes(x=Inf, y=Inf, #color="gray",
                                                      label=paste0(as.character(" \u03b2="),signif(beta,2),
                                                                   " [",signif(beta_min,2),", ",signif(beta_max,2),"],\nR\u00b2=",signif(R2,2))),
                                          hjust=1, vjust=1, show.legend=FALSE, 
                                          colour=mid.color
      )
      ggsave(paste0(d2,"_",d1,"_viability_log10_ratios_musyc.svg"),conf.plotB,width=9,height=4, device="svg")
      ggsave(paste0(d2,"_",d1,"_viability_log10_ratios_musyc.png"),conf.plotB,width=9,height=4, device="png")
      
      conf.plotB <- conf.plot + geom_text(data=combo.res, 
                                          mapping=aes(x=Inf, y=Inf, fontface=textFace,
                                                      label=paste0(as.character(" \u03b2="),signif(beta,2),
                                                                   " [",signif(beta_min,2),", ",signif(beta_max,2),"],\nR\u00b2=",signif(R2,2))),
                                          hjust=1, vjust=1, show.legend=FALSE, 
                                          colour=mid.color
      )
      ggsave(paste0(d2,"_",d1,"_viability_log10_ratios_allMusyc.svg"),conf.plotB,width=9,height=4, device="svg")
      ggsave(paste0(d2,"_",d1,"_viability_log10_ratios_allMusyc.png"),conf.plotB,width=9,height=4, device="png")
      
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
        conf.plotB <- conf.plot + geom_text(data=combo.resp, 
                                            mapping=aes(x=Inf, y=Inf, fontface=textFace,
                                                        label=paste0(as.character(" \u03b2="),signif(beta,2),
                                                                     " [",signif(beta_min,2),", ",signif(beta_max,2),"],\nR\u00b2=",signif(R2,2),
                                                                     "\np=",signif(p,2))),
                                            hjust=1, vjust=1, show.legend=FALSE, 
                                            colour=mid.color
        )
        ggsave(paste0(d2,"_",d1,"_viability_log10_ratios_allMusyc_p.svg"),conf.plotB,width=9,height=4, device="svg")
        ggsave(paste0(d2,"_",d1,"_viability_log10_ratios_allMusyc_p.png"),conf.plotB,width=9,height=4, device="png")
        
      }
    }
    #temp.colsA <- colorRampPalette(c(cols[[d2]], "gray"))(floor(length(unique(combo.res$conc12.ratio))/2))
    #temp.colsB <- colorRampPalette(c("gray",cols[[d1]]))(floor(length(unique(combo.res$conc12.ratio))/2))
    #temp.cols <- unique(c(temp.colsA, "gray", temp.colsB))
    temp.cols <- colorRampPalette(c(cols[[d2]], "gray", cols[[d1]]))(length(unique(combo.res$conc12.ratio)))
    conf.plot <- ggplot(combo.res, aes(x=combo.conc, y=100*meanGROWTH, color=as.factor(conc12.ratio))) +
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
    ggsave(paste0(d2,"_",d1,"_viability_log10_ratios_gray.pdf"),conf.plot,width=9,height=4)
    
    # conf.plotB <- conf.plot + geom_text(data=subset(combo.res, potentialSynergy), 
    #                                     mapping=aes(x=Inf, y=Inf, #color="gray", 
    #                                                 label=paste0("B=",signif(beta,2),
    #                                                              " [",signif(beta_min,2),", ",signif(beta_max,2),"],\nR2=",signif(R2,2))),
    #                                     hjust=1, vjust=1, show.legend=FALSE
    #                                     )
    #ggsave(paste0(d2,"_",d1,"_viability_log10_ratios_gray_musyc.pdf"),conf.plotB,width=9,height=4)
    
    if (any(combo.res$potentialSynergy)) {
      conf.plotB <- conf.plot + geom_text(data=subset(combo.res, potentialSynergy), 
                                          mapping=aes(x=Inf, y=Inf, #color="gray",
                                                      label=paste0(as.character(" \u03b2="),signif(beta,2),
                                                                   " [",signif(beta_min,2),", ",signif(beta_max,2),"],\nR\u00b2=",signif(R2,2))),
                                          hjust=1, vjust=1, show.legend=FALSE, 
                                          colour="gray")
      ggsave(paste0(d2,"_",d1,"_viability_log10_ratios_gray_musyc.svg"),conf.plotB,width=9,height=4, device="svg")
      ggsave(paste0(d2,"_",d1,"_viability_log10_ratios_gray_musyc.png"),conf.plotB,width=9,height=4, device="png")
    
      conf.plotB <- conf.plot + geom_text(data=combo.res, 
                                          mapping=aes(x=Inf, y=Inf, fontface=textFace,
                                                      label=paste0(as.character(" \u03b2="),signif(beta,2),
                                                                   " [",signif(beta_min,2),", ",signif(beta_max,2),"],\nR\u00b2=",signif(R2,2))),
                                          hjust=1, vjust=1, show.legend=FALSE, 
                                          colour="gray"
      )
      ggsave(paste0(d2,"_",d1,"_viability_log10_ratios_gray_allMusyc.svg"),conf.plotB,width=9,height=4, device="svg")
      ggsave(paste0(d2,"_",d1,"_viability_log10_ratios_gray_allMusyc.png"),conf.plotB,width=9,height=4, device="png")
    }
    
    temp.cols <- colorRampPalette(c(cols[[d1]], cols[[d2]]))(length(unique(combo.res$conc21.ratio)))
    conf.plot <- ggplot(combo.res, aes(x=combo.conc, y=100*meanGROWTH, color=as.factor(conc21.ratio))) +
      facet_grid(timeD ~ sample) +
      geom_smooth(linetype="dashed", se=FALSE, method=drc::drm, method.args=list(fct=L.4(), se=FALSE)) +
      #scale_color_manual(values=c(scales::hue_pal()(2))) + 
      #scale_shape_manual(values=1:length(all.mpnst), breaks=all.mpnst)+
      scale_x_continuous(transform="log10") + geom_point() + 
      geom_errorbar(aes(ymin=100*(meanGROWTH-sdGROWTH), ymax=100*(meanGROWTH+sdGROWTH))) +
      #ggrepel::geom_text_repel(combo.res, aes(x=0,y=0,label=title))+
      ggtitle(paste0(d1," + ",d2))+scale_color_manual(values=temp.cols) +
      theme_classic(base_size=12) + labs(x="Concentration (uM)", y = "% Viability",
                                         color = paste("Ratio of",d2,"\nto",d1)) +
      theme(plot.title=element_text(hjust=0.5,face="bold"), 
            axis.text.x=element_text(angle=45, hjust=1, vjust=1))
    ggsave(paste0(d2,"_",d1,"_viability_log10_ratios_v2.pdf"),conf.plot,width=9,height=4)
    
    # conf.plotB <- conf.plot + geom_text(data=subset(combo.res, potentialSynergy), 
    #                                     mapping=aes(x=Inf, y=Inf, #color="gray",
    #                                                 label=paste0("B=",signif(beta,2),
    #                                                              " [",signif(beta_min,2),", ",signif(beta_max,2),"],\nR2=",signif(R2,2))),
    #                                     hjust=1, vjust=1, show.legend=FALSE
    #                                     )
    # ggsave(paste0(d2,"_",d1,"_viability_log10_ratios_v2_musyc.pdf"),conf.plotB,width=9,height=4)
    # 
    if (any(combo.res$potentialSynergy)) {
      conf.plotB <- conf.plot + geom_text(data=subset(combo.res, potentialSynergy), 
                                          mapping=aes(x=Inf, y=Inf, #color="gray",
                                                      label=paste0(as.character(" \u03b2="),signif(beta,2),
                                                                   " [",signif(beta_min,2),", ",signif(beta_max,2),"],\nR\u00b2=",signif(R2,2))),
                                          hjust=1, vjust=1, show.legend=FALSE, 
                                          colour=mid.color)
      ggsave(paste0(d2,"_",d1,"_viability_log10_ratios_v2_musyc.svg"),conf.plotB,width=9,height=4, device="svg")
      ggsave(paste0(d2,"_",d1,"_viability_log10_ratios_v2_musyc.png"),conf.plotB,width=9,height=4, device="png")
      
      conf.plotB <- conf.plot + geom_text(data=combo.res, 
                                          mapping=aes(x=Inf, y=Inf, fontface=textFace,
                                                      label=paste0(as.character(" \u03b2="),signif(beta,2),
                                                                   " [",signif(beta_min,2),", ",signif(beta_max,2),"],\nR\u00b2=",signif(R2,2))),
                                          hjust=1, vjust=1, show.legend=FALSE, 
                                          colour=mid.color
      )
      ggsave(paste0(d2,"_",d1,"_viability_log10_ratios_v2_allMusyc.svg"),conf.plotB,width=9,height=4, device="svg")
      ggsave(paste0(d2,"_",d1,"_viability_log10_ratios_v2_allMusyc.png"),conf.plotB,width=9,height=4, device="png")
    }
    
    #temp.colsA <- colorRampPalette(c(cols[[d1]], "gray"))(floor(length(unique(combo.res$conc21.ratio))/2))
    #temp.colsB <- colorRampPalette(c("gray",cols[[d2]]))(floor(length(unique(combo.res$conc21.ratio))/2))
    #temp.cols <- unique(c(temp.colsA, "gray", temp.colsB))
    temp.cols <- colorRampPalette(c(cols[[d1]], "gray", cols[[d2]]))(length(unique(combo.res$conc21.ratio)))
    conf.plot <- ggplot(combo.res, aes(x=combo.conc, y=100*meanGROWTH, color=as.factor(conc21.ratio))) +
      facet_grid(timeD ~ sample) +
      geom_smooth(linetype="dashed", se=FALSE, method=drc::drm, method.args=list(fct=L.4(), se=FALSE)) +
      #scale_color_manual(values=c(scales::hue_pal()(2))) + 
      #scale_shape_manual(values=1:length(all.mpnst), breaks=all.mpnst)+
      scale_x_continuous(transform="log10") + geom_point() + 
      geom_errorbar(aes(ymin=100*(meanGROWTH-sdGROWTH), ymax=100*(meanGROWTH+sdGROWTH))) +
      #ggrepel::geom_text_repel(combo.res, aes(x=0,y=0,label=title))+
      ggtitle(paste0(d1," + ",d2))+scale_color_manual(values=temp.cols) +
      theme_classic(base_size=12) + labs(x="Concentration (uM)", y = "% Viability",
                                         color = paste("Ratio of",d2,"\nto",d1)) +
      theme(plot.title=element_text(hjust=0.5,face="bold"), 
            axis.text.x=element_text(angle=45, hjust=1, vjust=1))
    ggsave(paste0(d2,"_",d1,"_viability_log10_ratios_v2_gray.pdf"),conf.plot,width=9,height=4)

    # conf.plotB <- conf.plot + geom_text(data=subset(combo.res, potentialSynergy), 
    #                                     mapping=aes(x=Inf, y=Inf, #color="gray",
    #                                                 label=paste0("B=",signif(beta,2),
    #                                                              " [",signif(beta_min,2),", ",signif(beta_max,2),"],\nR2=",signif(R2,2))),
    #                                     hjust=1, vjust=1, show.legend=FALSE
    #                                     )
    # ggsave(paste0(d2,"_",d1,"_viability_log10_ratios_v2_gray_musyc.pdf"),conf.plotB,width=9,height=4)
    # 
    if (any(combo.res$potentialSynergy)) {
      conf.plotB <- conf.plot + geom_text(data=subset(combo.res, potentialSynergy), 
                                          mapping=aes(x=Inf, y=Inf, #color="gray",
                                                      label=paste0(as.character(" \u03b2="),signif(beta,2),
                                                                   " [",signif(beta_min,2),", ",signif(beta_max,2),"],\nR\u00b2=",signif(R2,2))),
                                          hjust=1, vjust=1, show.legend=FALSE, 
                                          colour="gray")
      ggsave(paste0(d2,"_",d1,"_viability_log10_ratios_v2_gray_musyc.svg"),conf.plotB,width=9,height=4, device="svg")
      ggsave(paste0(d2,"_",d1,"_viability_log10_ratios_v2_gray_musyc.png"),conf.plotB,width=9,height=4, device="png")
      
      conf.plotB <- conf.plot + geom_text(data=combo.res, 
                                          mapping=aes(x=Inf, y=Inf, fontface=textFace,
                                                      label=paste0(as.character(" \u03b2="),signif(beta,2),
                                                                   " [",signif(beta_min,2),", ",signif(beta_max,2),"],\nR\u00b2=",signif(R2,2))),
                                          hjust=1, vjust=1, show.legend=FALSE, 
                                          colour="gray"
      )
      ggsave(paste0(d2,"_",d1,"_viability_log10_ratios_v2_gray_allMusyc.svg"),conf.plotB,width=9,height=4, device="svg")
      ggsave(paste0(d2,"_",d1,"_viability_log10_ratios_v2_gray_allMusyc.png"),conf.plotB,width=9,height=4, device="png")
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
shared.corr <- na.omit(shared.corr[,c("sample","drugCombo",
                                      "Pearson.est","Pearson.p","Pearson.q",
                                      "Spearman.est","Spearman.p","Spearman.q","N")])
drug.corr.res <- merge(shared.corr, results, by=c("drugCombo","sample")) # no time overlap

# note: only mirdametinib+doxorubicin is significantly correlated in MN-2 (Pearson.q<0.05) which has Pearson.est=0.1915
musyc.drug.corr <- cor.test(drug.corr.res$beta, drug.corr.res$Pearson.est) # r: 0.2312464, p = 2.760436e-14
musyc.drug.corr$estimate
musyc.drug.corr$p.value
musyc.drug.corrsp <- cor.test(drug.corr.res$beta, drug.corr.res$Pearson.est,method="spearman") # rho: 0.07683617, p = 0.01250316
musyc.drug.corrsp$estimate
musyc.drug.corrsp$p.value

musyc.drug.corr2 <- cor.test(drug.corr.res$beta, drug.corr.res$Spearman.est) # r: 0.1911145, p = 3.827007e-10
musyc.drug.corr2$estimate
musyc.drug.corr2$p.value
musyc.drug.corr2sp <- cor.test(drug.corr.res$beta, drug.corr.res$Spearman.est,method="spearman") # rho: 0.07740114, p = 0.01186832
musyc.drug.corr2sp$estimate
musyc.drug.corr2sp$p.value

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

ggplot(drug.corr.res, aes(x=beta, y=Pearson.est, #shape=sample, # only sample is JH-2-002
                          color=color, label=drugCombo)) + 
  geom_point() + theme_minimal() + scale_color_identity() + ggrepel::geom_label_repel() +
  labs(x="MuSyC: % Change in Viability with Combination vs. Most Effective Single Agent",
       y="DMEA: Predicted Drug Sensitivity")
shared.combos2 <- unique(drug.corr.res$drugCombo) # 11

drug.corr.res2 <- na.omit(drug.corr.res[,c("drugCombo","sample","beta",
                                           "Pearson.est","Pearson.p","Pearson.q",
                                           "Spearman.est","Spearman.p","Spearman.q","color")])
ggplot(drug.corr.res2, aes(x=beta, y=Pearson.est, #shape=sample, # only sample is JH-2-002
                          color=drugCombo)) + 
  geom_point() + theme_minimal() + scale_color_manual(values=combo.cols, breaks=names(combo.cols), labels=names(combo.cols)) + 
  labs(x="MuSyC: % Change in Viability with Combination vs. Most Effective Single Agent",
       y="DMEA: Predicted Drug Sensitivity")
ggplot(drug.corr.res2, aes(x=beta, y=Spearman.est, #shape=sample, # only sample is JH-2-002
                           color=drugCombo)) + 
  geom_point() + theme_minimal() + scale_color_manual(values=combo.cols, breaks=names(combo.cols), labels=names(combo.cols)) + 
  labs(x="MuSyC: % Change in Viability with Combination vs. Most Effective Single Agent",
       y="DMEA: Predicted Drug Sensitivity")

ggplot(drug.corr.res2[abs(drug.corr.res2$beta)<2,], aes(x=beta, y=Pearson.est, #shape=sample, # only sample is JH-2-002
                           color=drugCombo)) + 
  geom_point() + theme_minimal() + scale_color_manual(values=combo.cols, breaks=names(combo.cols), labels=names(combo.cols)) + 
  labs(x="MuSyC: % Change in Viability with Combination vs. Most Effective Single Agent",
       y="DMEA: Predicted Drug Sensitivity")
musyc.drug.corr <- cor.test(drug.corr.res2[abs(drug.corr.res2$beta)<2,]$beta, drug.corr.res2[abs(drug.corr.res2$beta)<2,]$Pearson.est)
musyc.drug.corr$estimate # 0.04058664
musyc.drug.corr$p.value # 0.2089673
musyc.drug.corrsp <- cor.test(drug.corr.res2[abs(drug.corr.res2$beta)<2,]$beta, drug.corr.res2[abs(drug.corr.res2$beta)<2,]$Pearson.est,method="spearman")
musyc.drug.corrsp$estimate # 0.02069031 
musyc.drug.corrsp$p.value # 0.5219781

ggplot(drug.corr.res2[abs(drug.corr.res2$beta)<2,], aes(x=beta, y=Spearman.est, #shape=sample, # only sample is JH-2-002
                                                    color=drugCombo)) + 
  geom_point() + theme_minimal() + scale_color_manual(values=combo.cols, breaks=names(combo.cols), labels=names(combo.cols)) + 
  labs(x="MuSyC: % Change in Viability with Combination vs. Most Effective Single Agent",
       y="DMEA: Predicted Drug Sensitivity")
musyc.drug.corr <- cor.test(drug.corr.res2[abs(drug.corr.res2$beta)<2,]$beta, drug.corr.res2[abs(drug.corr.res2$beta)<2,]$Spearman.est)
musyc.drug.corr$estimate # 0.02122183
musyc.drug.corr$p.value # 0.5113407
musyc.drug.corrsp <- cor.test(drug.corr.res2[abs(drug.corr.res2$beta)<2,]$beta, drug.corr.res2[abs(drug.corr.res2$beta)<2,]$Spearman.est,method="spearman")
musyc.drug.corrsp$estimate # 0.0150475  
musyc.drug.corrsp$p.value # 0.641467

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

results2 <- results[,c("sample","drug1","drug2","beta","R2","timeD")]
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

# ggplot(shared.moa.res, aes(x=beta, y=NES, #shape=sample, # only sample is JH-2-002
#                                                         color=moaComboShort)) + 
#   geom_point() + theme_minimal() + #scale_color_manual(values=combo.cols, breaks=names(combo.cols), labels=names(combo.cols)) + 
#   labs(x="MuSyC: % Change in Viability with Combination vs. Most Effective Single Agent",
#        y="DMEA: Predicted Drug Sensitivity")
# musyc.drug.corr <- cor.test(shared.moa.res$beta, shared.moa.res$NES)
# musyc.drug.corr$estimate # 0.2290539
# musyc.drug.corr$p.value # 2.706572e-08
# musyc.drug.corrsp <- cor.test(shared.moa.res$beta, shared.moa.res$NES ,method="spearman")
# musyc.drug.corrsp$estimate # 0.2332287 
# musyc.drug.corrsp$p.value # 1.482147e-08
# 
# ggplot(shared.moa.res[shared.moa.res$sig,], aes(x=beta, y=NES, #shape=sample, # only sample is JH-2-002
#                            color=moaComboShort)) + 
#   geom_point() + theme_minimal() + #scale_color_manual(values=combo.cols, breaks=names(combo.cols), labels=names(combo.cols)) + 
#   labs(x="MuSyC: % Change in Viability with Combination vs. Most Effective Single Agent",
#        y="DMEA: Predicted Drug Sensitivity")
# musyc.drug.corr <- cor.test(shared.moa.res[shared.moa.res$sig,]$beta, shared.moa.res[shared.moa.res$sig,]$NES)
# musyc.drug.corr$estimate # 0.1979562 
# musyc.drug.corr$p.value # 0.000729327
# musyc.drug.corrsp <- cor.test(shared.moa.res[shared.moa.res$sig,]$beta, shared.moa.res[shared.moa.res$sig,]$NES ,method="spearman")
# musyc.drug.corrsp$estimate # 0.1309307 
# musyc.drug.corrsp$p.value # 0.02629154
# 
# ggplot(shared.moa.res[abs(shared.moa.res$beta)<2,], aes(x=beta, y=NES)) + 
#   geom_point(aes(shape=moaComboShort, color=sig, size=R2)) + 
#   theme_minimal() + ggrepel::geom_text_repel(aes(label=paste(as.character(Time), "RNA,",as.character(timeD),"viability"))) + 
#   #coord_flip()+stat_smooth(method="lm", se=FALSE, formula=y ~ x + I(x^2)) +
#   labs(x="MuSyC: % Change in Viability with Combination vs. Most Effective Single Agent",
#        y="DMEA: Predicted Drug Sensitivity")
# 
# ggplot(shared.moa.res[abs(shared.moa.res$beta)<2,], aes(x=beta, y=NES)) + 
#   geom_point(aes(shape=moaComboShort, color=sig#, size=R2)
#                  )) + 
#   theme_minimal() + ggrepel::geom_text_repel(aes(label=paste(as.character(Time), "RNA,",as.character(timeD),"viability"))) + 
#   geom_smooth(method="lm",formula = x ~ poly(y,2), orientation="y")+
#   #coord_flip()+stat_smooth(method="lm", se=FALSE, formula=y ~ x + I(x^2)) +
#   labs(x="MuSyC: % Change in Viability with Combination vs. Most Effective Single Agent",
#        y="DMEA: Predicted Drug Sensitivity")
# 
shared.moa.res2 <- dplyr::distinct(shared.moa.res)
# ggplot(shared.moa.res2[abs(shared.moa.res2$beta)<2,], aes(y=beta, x=NES)) + 
#   geom_point(aes(shape=moaComboShort, color=sig)) + 
#   theme_minimal() + ggrepel::geom_text_repel(data=subset(shared.moa.res2[abs(shared.moa.res2$beta)<2,],sig),
#                                              aes(label=paste(as.character(Time), "RNA,",as.character(timeD),"viability"))) + 
#   geom_smooth(method="lm",formula = y ~ poly(x,2))+
#   labs(y="MuSyC: % Change in Viability with Combination vs. Most Effective Single Agent",
#        x="DMEA: Predicted Drug Sensitivity")
# 
# ggplot(shared.moa.res2[abs(shared.moa.res2$beta)<2 & shared.moa.res2$beta<0,], aes(y=beta, x=NES)) + 
#   geom_point(aes(shape=moaComboShort, color=sig)) + 
#   theme_minimal() + ggrepel::geom_text_repel(data=subset(shared.moa.res2[abs(shared.moa.res2$beta)<2 &
#                                                                            shared.moa.res2$beta<0,],sig),
#                                              aes(label=paste(as.character(Time), "RNA,",as.character(timeD),"viability"))) + 
#   geom_smooth(method="lm",formula = y ~ poly(x,2))+
#   labs(y="MuSyC: % Change in Viability with Combination vs. Most Effective Single Agent",
#        x="DMEA: Predicted Drug Sensitivity")
# 
# ggplot(shared.moa.res2, aes(y=beta, x=NES)) + 
#   geom_point(aes(shape=moaComboShort, color=sig)) + 
#   theme_minimal() + ggrepel::geom_text_repel(data=subset(shared.moa.res2,sig),
#                                              aes(label=paste(as.character(Time), "RNA,",as.character(timeD),"viability"))) + 
#   geom_smooth(method="lm",formula = y ~ poly(x,2))+
#   labs(y="MuSyC: % Change in Viability with Combination vs. Most Effective Single Agent",
#        x="DMEA: Predicted Drug Sensitivity")
# 
# ggplot(shared.moa.res2, aes(y=-beta, x=NES)) + 
#   geom_point(aes(shape=moaComboShort, color=sig)) + 
#   theme_minimal() + ggrepel::geom_text_repel(data=subset(shared.moa.res2,sig),
#                                              aes(label=paste(as.character(Time), "RNA,",as.character(timeD),"viability"))) + 
#   scale_y_continuous(transform="log10") +
#   #geom_smooth(method="lm",formula = y ~ poly(x,2))+ #scale_y_continuous(transform="log10") +
#   labs(y="MuSyC: -% Change in Viability with Combination",
#        x="DMEA: Predicted Drug Sensitivity", color = "Significant Enrichment") + 
#   scale_color_manual(values=c("red","grey"), breaks=c(TRUE, FALSE))
# 
# ggplot(shared.moa.res2, aes(y=-beta, x=NES)) + 
#   geom_point(aes(shape=moaComboShort, color=sig)) + 
#   theme_minimal() + 
#   scale_y_continuous(transform="log10") +
#   geom_smooth(method="lm",se=FALSE,formula = y ~ poly(x,2))+ #scale_y_continuous(transform="log10") +
#   labs(y="MuSyC Synergy Score", shape = "Drug Combination",
#        x="DMEA Sensitivity Score", color = "Significant\nDMEA Result") + 
#   scale_color_manual(values=c("red","grey"), breaks=c(TRUE, FALSE))
# 
# ggplot(shared.moa.res2, aes(y=-beta, x=NES)) + 
#   geom_point(aes(shape=moaComboShort, color=sig)) + 
#   theme_minimal() + 
#   scale_y_continuous(transform="log10") +
#   geom_smooth(method="lm",se=FALSE,formula = y ~ x + I(x^2), size=1)+ 
#   labs(y="MuSyC Synergy Score", shape = "Drug Combination",
#        x="DMEA Sensitivity Score", color = "Significant\nDMEA Result") + 
#   scale_color_manual(values=c("red","grey"), breaks=c(TRUE, FALSE))
# 
# ggplot(shared.moa.res2, aes(y=-beta, x=NES)) + 
#   geom_point(aes(shape=moaComboShort, color=sig)) + 
#   theme_minimal() + 
#   scale_y_continuous(transform="log2") +#scale_y_continuous(transform="-") +
#   geom_smooth(method="lm",se=FALSE,formula = y ~ x + I(x^2), size=1)+ 
#   labs(y="MuSyC Synergy Score", shape = "Drug Combination",
#        x="DMEA Sensitivity Score", color = "Significant\nDMEA Result") + 
#   scale_color_manual(values=c("red","grey"), breaks=c(TRUE, FALSE))
# 
# ggplot(shared.moa.res2, aes(y=-beta, x=NES)) + 
#   geom_point(aes(shape=moaComboShort, color=sig)) + 
#   theme_minimal() + 
#   scale_y_continuous(transform="log2") +#scale_y_continuous(transform="-") +
#   #geom_smooth(method="lm",se=FALSE,formula = y ~ x + I(x^2), size=1)+ 
#   labs(y="MuSyC Synergy Score", shape = "Drug Combination",
#        x="DMEA Sensitivity Score", color = "Significant\nDMEA Result") + 
#   scale_color_manual(values=c("red","grey"), breaks=c(TRUE, FALSE))
# 
# ggplot(shared.moa.res2, aes(y=-beta, x=NES)) + 
#   geom_point(aes(shape=moaComboShort, color=sig, size=R2)) + 
#   theme_minimal() + ggrepel::geom_text_repel(data=subset(shared.moa.res2,sig),
#                                              aes(label=paste(as.character(Time), "RNA,\n",as.character(timeD),"viability"))) + 
#   labs(y="MuSyC Synergy Score", shape = "Drug Combination", size=expression("R"^2),
#        x="DMEA Sensitivity Score", color = "Significant\nDMEA Result") + 
#   scale_color_manual(values=c("red","grey"), breaks=c(TRUE, FALSE)) + scale_size_continuous(range=c(2,4))

ggplot(shared.moa.res2, aes(y=-beta, x=NES)) + 
  geom_point(aes(shape=moaComboShort, color=sig), size=3) + 
  theme_minimal() + ggrepel::geom_text_repel(data=subset(shared.moa.res2,sig),
                                             aes(label=paste(as.character(Time), "RNA,\n",as.character(timeD),"viability"))) + 
  labs(y="MuSyC Synergy Score", shape = "Drug Combination", 
       x="DMEA Sensitivity Score", color = "Significant\nDMEA Result") + 
  scale_color_manual(values=c("red","grey"), breaks=c(TRUE, FALSE))
ggsave("MOA_results_musyc_vs_dmea.pdf", width=7, height=4)

ggplot(shared.moa.res2, aes(y=-beta, x=NES)) + 
  geom_point(aes(shape=moaComboShort, color=sig), size=3) + 
  geom_smooth(method="lm",se=FALSE,formula = y ~ x + I(x^2), color="darkgrey", linetype="dashed") +
  theme_minimal() + ggrepel::geom_text_repel(data=subset(shared.moa.res2,sig),
                                             aes(label=paste(as.character(Time), "RNA,\n",as.character(timeD),"viability"))) + 
  labs(y="MuSyC Synergy Score", shape = "Drug Combination", 
       x="DMEA Sensitivity Score", color = "Significant\nDMEA Result") +
  scale_color_manual(values=c("red","grey"), breaks=c(TRUE, FALSE))
ggsave("MOA_results_musyc_vs_dmea_wQuadraticFit.pdf", width=7, height=4)

ggplot(shared.moa.res2, aes(y=-beta, x=NES)) + 
  geom_point(aes(shape=moaComboShort, color=sig), size=3) + 
  theme_minimal() + ggrepel::geom_text_repel(data=subset(shared.moa.res2,sig),
                                             aes(label=paste(as.character(Time), "RNA,\n",as.character(timeD),"viability"))) + 
  labs(y="MuSyC Synergy Score", shape = "Drug Combination", 
       x="DMEA Sensitivity Score", color = "Significant\nDMEA Result") + 
  scale_color_manual(values=c("red","grey"), breaks=c(TRUE, FALSE)) + scale_y_continuous(transform="log10")
ggsave("MOA_results_musyc_vs_dmea_log10.pdf", width=7, height=4)
write.csv(shared.moa.res2,"MOA_results_musyc_vs_dmea.csv", row.names=FALSE)

ggplot(shared.moa.res2[abs(shared.moa.res2$beta)<2,], aes(y=-beta, x=NES)) + 
  geom_point(aes(shape=moaComboShort, color=sig), size=3) + 
  theme_minimal() + ggrepel::geom_text_repel(data=subset(shared.moa.res2[abs(shared.moa.res2$beta)<2,],sig),
                                             aes(label=paste(as.character(Time), "RNA,\n",as.character(timeD),"viability"))) + 
  labs(y="MuSyC Synergy Score", shape = "Drug Combination", 
       x="DMEA Sensitivity Score", color = "Significant\nDMEA Result") + 
  scale_color_manual(values=c("red","grey"), breaks=c(TRUE, FALSE))
ggsave("MOA_results_musyc_vs_dmea_absMax2B.pdf", width=7, height=4)

ggplot(shared.moa.res2[abs(shared.moa.res2$beta)<2,], aes(y=-beta, x=NES)) + 
  geom_point(aes(shape=moaComboShort, color=sig), size=3) + 
  geom_smooth(method="lm",se=FALSE,formula = y ~ x + I(x^2), color="darkgrey", linetype="dashed") +
  theme_minimal() + ggrepel::geom_text_repel(data=subset(shared.moa.res2[abs(shared.moa.res2$beta)<2,],sig),
                                             aes(label=paste(as.character(Time), "RNA,\n",as.character(timeD),"viability"))) + 
  labs(y="MuSyC Synergy Score", shape = "Drug Combination", 
       x="DMEA Sensitivity Score", color = "Significant\nDMEA Result") +
  scale_color_manual(values=c("red","grey"), breaks=c(TRUE, FALSE))
ggsave("MOA_results_musyc_vs_dmea_wQuadraticFit_absMax2B.pdf", width=7, height=4)

ggplot(shared.moa.res2[abs(shared.moa.res2$beta)<2,], aes(y=-beta, x=NES)) + 
  geom_point(aes(shape=moaComboShort, color=sig), size=3) + 
  theme_minimal() + ggrepel::geom_text_repel(data=subset(shared.moa.res2[abs(shared.moa.res2$beta)<2,],sig),
                                             aes(label=paste(as.character(Time), "RNA,\n",as.character(timeD),"viability"))) + 
  labs(y="MuSyC Synergy Score", shape = "Drug Combination", 
       x="DMEA Sensitivity Score", color = "Significant\nDMEA Result") + 
  scale_color_manual(values=c("red","grey"), breaks=c(TRUE, FALSE)) + scale_y_continuous(transform="log10")
ggsave("MOA_results_musyc_vs_dmea_log10_absMax2B.pdf", width=7, height=4)

ggplot(shared.moa.res2[abs(shared.moa.res2$beta)<2 & shared.moa.res2$beta<0,], aes(y=-beta, x=NES)) + 
  geom_point(aes(shape=moaComboShort, color=sig), size=3) + 
  theme_minimal() + ggrepel::geom_text_repel(data=subset(shared.moa.res2[abs(shared.moa.res2$beta)<2 & shared.moa.res2$beta<0,],sig),
                                             aes(label=paste(as.character(Time), "RNA,\n",as.character(timeD),"viability"))) + 
  labs(y="MuSyC Synergy Score", shape = "Drug Combination", 
       x="DMEA Sensitivity Score", color = "Significant\nDMEA Result") + 
  scale_color_manual(values=c("red","grey"), breaks=c(TRUE, FALSE))
ggsave("MOA_results_musyc_vs_dmea_absMax2B_max0B.pdf", width=7, height=4)

ggplot(shared.moa.res2[abs(shared.moa.res2$beta)<2 & shared.moa.res2$beta<0,], aes(y=-beta, x=NES)) + 
  geom_point(aes(shape=moaComboShort, color=sig), size=3) + 
  geom_smooth(method="lm",se=FALSE,formula = y ~ x + I(x^2), color="darkgrey", linetype="dashed") +
  theme_minimal() + ggrepel::geom_text_repel(data=subset(shared.moa.res2[abs(shared.moa.res2$beta)<2 & shared.moa.res2$beta<0,],sig),
                                             aes(label=paste(as.character(Time), "RNA,\n",as.character(timeD),"viability"))) + 
  labs(y="MuSyC Synergy Score", shape = "Drug Combination", 
       x="DMEA Sensitivity Score", color = "Significant\nDMEA Result") +
  scale_color_manual(values=c("red","grey"), breaks=c(TRUE, FALSE))
ggsave("MOA_results_musyc_vs_dmea_wQuadraticFit_absMax2B_max0B.pdf", width=7, height=4)

ggplot(shared.moa.res2[abs(shared.moa.res2$beta)<2 & shared.moa.res2$beta<0,], aes(y=-beta, x=NES)) + 
  geom_point(aes(shape=moaComboShort, color=sig), size=3) + 
  theme_minimal() + ggrepel::geom_text_repel(data=subset(shared.moa.res2[abs(shared.moa.res2$beta)<2 & shared.moa.res2$beta<0,],sig),
                                             aes(label=paste(as.character(Time), "RNA,\n",as.character(timeD),"viability"))) + 
  labs(y="MuSyC Synergy Score", shape = "Drug Combination", 
       x="DMEA Sensitivity Score", color = "Significant\nDMEA Result") + 
  scale_color_manual(values=c("red","grey"), breaks=c(TRUE, FALSE)) + scale_y_continuous(transform="log10")
ggsave("MOA_results_musyc_vs_dmea_log10_absMax2B_max0B.pdf", width=7, height=4)

ggplot(shared.moa.res2, aes(y=-beta, x=NES)) + 
  geom_point(aes(shape=moaComboShort, color=sig), size=3) + 
  theme_minimal() + ggrepel::geom_text_repel(data=subset(shared.moa.res2,sig),
                                             aes(label=paste(as.character(Time), "RNA,\n",as.character(timeD),"viability"))) + 
  labs(y="MuSyC Synergy Score", shape = "Drug Combination", 
       x="DMEA Sensitivity Score", color = "Significant\nDMEA Result") + 
  scale_color_manual(values=c("red","grey"), breaks=c(TRUE, FALSE)) + scale_y_continuous(transform="log10") +
  geom_smooth(method="lm",se=FALSE,formula = y ~ x + I(x^2), color="darkgrey", linetype="dashed")
ggsave("MOA_results_musyc_vs_dmea_log10_wQuadtraticFit.pdf", width=7, height=4)

ggplot(shared.moa.res2, aes(y=-beta, x=NES)) + 
  geom_point(aes(shape=moaComboShort, color=sig), size=3) + 
  theme_minimal() + ggrepel::geom_text_repel(data=subset(shared.moa.res2,sig),
                                             aes(label=paste(as.character(Time), "RNA,\n",as.character(timeD),"viability"))) + 
  labs(y="MuSyC Synergy Score", shape = "Drug Combination", 
       x="DMEA Sensitivity Score", color = "Significant\nDMEA Result") + 
  scale_color_manual(values=c("red","grey"), breaks=c(TRUE, FALSE)) + scale_y_continuous(transform="log2") +
  geom_smooth(method="lm",se=FALSE,formula = y ~ x + I(x^2), color="darkgrey", linetype="dashed")
ggsave("MOA_results_musyc_vs_dmea_log2_wQuadtraticFit.pdf", width=7, height=4)

# are beta values lower in significant DMEA results?
b.test <- t.test(shared.moa.res2[shared.moa.res2$sig,]$beta, shared.moa.res2[!shared.moa.res2$sig,]$beta, alternative="less")
b.test$p.value # 0.1032556
b.test$estimate # mean of x mean of y 
# -2.617753 -1.026053 
length(shared.moa.res2[shared.moa.res2$sig,]$beta) # 12
length(shared.moa.res2[!shared.moa.res2$sig,]$beta) # 12
# 
# ggplot(shared.moa.res[abs(shared.moa.res$beta)<2,], aes(x=beta, y=NES)) + 
#   geom_point(aes(shape=moaComboShort, color=sig)) + theme_minimal() + 
#   ggrepel::geom_text_repel(aes(label=paste(as.character(Time), "RNA,",as.character(timeD),"viability"))) + 
#   labs(x="MuSyC: % Change in Viability with Combination vs. Most Effective Single Agent",
#        y="DMEA: Predicted Drug Sensitivity")+
#   coord_flip()#+#stat_smooth(method="lm", se=FALSE, formula=y ~ x + I(x^2))
#   
# ggplot(shared.moa.res[abs(shared.moa.res$beta)<2,], aes(y=beta, x=NES)) + 
#   geom_point(aes(shape=moaComboShort, color=sig, size=R2)) + 
#   theme_minimal() + ggrepel::geom_text_repel(aes(label=paste(as.character(Time), "RNA,",as.character(timeD),"viability"))) + 
#   stat_smooth(method="lm", se=FALSE, formula=y ~ x^2) +
#   labs(y="MuSyC: % Change in Viability with Combination vs. Most Effective Single Agent",
#        x="DMEA: Predicted Drug Sensitivity")
# musyc.drug.corr <- cor.test(shared.moa.res[abs(shared.moa.res$beta)<2,]$beta, shared.moa.res[abs(shared.moa.res$beta)<2,]$NES)
# musyc.drug.corr$estimate # 0.1519407
# musyc.drug.corr$p.value # 0.0008390425
# musyc.drug.corrsp <- cor.test(shared.moa.res[abs(shared.moa.res$beta)<2,]$beta, shared.moa.res[abs(shared.moa.res$beta)<2,]$NES,method="spearman")
# musyc.drug.corrsp$estimate # 0.1168484  
# musyc.drug.corrsp$p.value # 0.01040329
# 
# ggplot(shared.moa.res[abs(shared.moa.res$beta)<2 & shared.moa.res$sig,], aes(x=beta, y=NES, #shape=sample, # only sample is JH-2-002
#                                                 color=moaComboShort)) + 
#   geom_point() + theme_minimal() + 
#   labs(x="MuSyC: % Change in Viability with Combination vs. Most Effective Single Agent",
#        y="DMEA: Predicted Drug Sensitivity")
# musyc.drug.corr <- cor.test(shared.moa.res[abs(shared.moa.res$beta)<2 & shared.moa.res$sig,]$beta, 
#                             shared.moa.res[abs(shared.moa.res$beta)<2 & shared.moa.res$sig,]$NES)
# musyc.drug.corr$estimate # 0.1322933 
# musyc.drug.corr$p.value # 0.05218954
# musyc.drug.corrsp <- cor.test(shared.moa.res[abs(shared.moa.res$beta)<2 & shared.moa.res$sig,]$beta, 
#                               shared.moa.res[abs(shared.moa.res$beta)<2 & shared.moa.res$sig,]$NES ,method="spearman")
# musyc.drug.corrsp$estimate # 0 
# musyc.drug.corrsp$p.value # 1
