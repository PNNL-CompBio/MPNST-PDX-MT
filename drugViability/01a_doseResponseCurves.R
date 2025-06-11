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
#### drug1 centric ####
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

#### drug2 centric ####
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

#### ratio centric ####
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
