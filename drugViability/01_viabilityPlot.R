# viability plot
library(plyr);library(dplyr);library(synapser);library(ggplot2)
synapser::synLogin()
base.path <- "~/OneDrive - PNNL/Documents/GitHub/MPNST-PDX-MT/drugViability"
setwd(base.path)

# run 00_createCurveStats.py
library(reticulate)
source_python("00_createCurveStats.py")

drug.info <- list("palbociclib" = "CDK inhibitor", "ribociclib" = "CDK inhibitor",
                  "trabectedin" = "Chemotherapy", "ifosfamide" = "DNA alkylating agent",
                  "decitabine" = "DNMT inhibitor", "vorinostat" = "HDAC inhibitor",
                  "mirdametinib" = "MEK inhibitor", "selumetinib" = "MEK inhibitor",
                  "trametinib" = "MEK inhibitor", "capmatinib" = "MET inhibitor",
                  "olaparib" = "PARP inhibitor", "RMC4630" = "SHP2 inhibitor",
                  "TNO155" = "SHP2 inhibitor", "doxorubicin" = "TOP inhibitor",
                  "irinotecan" = "TOP inhibitor", "verteporfin" = "YAP inhibitor")

## load data
viability <- read.table(synapser::synGet("syn65941820")$path, sep="\t", header = TRUE)
viability <- viability[viability$improve_drug_id != "irinotecan",]
viability.cmax <- readxl::read_excel(synapser::synGet("syn70107593")$path)

## plot
# viability at Cmax
mean.viability.cmax <- plyr::ddply(viability.cmax, .(Drug, PDX, Hour), summarize,
                                   meanViability = mean(`Viability at Cmax`, na.rm=TRUE),
                                   sdViability = sd(`Viability at Cmax`, na.rm=TRUE))
mean.viability.cmax$timeD <- paste0(as.numeric(mean.viability.cmax$Hour)/24,"d")
# need to scale SD such that bigger SD is associated with smaller point sizes
mean.viability.cmax$size <- 1/mean.viability.cmax$sdViability
ggplot(mean.viability.cmax, aes(x=PDX, y=reorder(Drug, -meanViability),
                   color=meanViability, size=size)) + 
  geom_point() + 
  facet_wrap(.~timeD)+theme_classic() + 
  scale_color_gradient(high="grey",low="red",limits=c(0,ceiling(max(mean.viability.cmax$meanViability)))) + 
  theme(axis.title = element_blank(), axis.text.x = element_text(angle=45, vjust=1, hjust=1)) + 
  labs(size="1/SD",color="Mean Viability\nat Cmax")
ggsave("CmaxViability_heatmap_orderMeanViability_sizeInverseSD.pdf", width=4,height=4)

ggplot(mean.viability.cmax, aes(x=PDX, y=reorder(Drug, -meanViability),
                                color=meanViability, size=sdViability^-0.5)) + 
  geom_point() + 
  facet_wrap(.~timeD)+theme_classic() + 
  scale_color_gradient(high="grey",low="red",limits=c(0,ceiling(max(mean.viability.cmax$meanViability)))) + 
  theme(axis.title = element_blank(), axis.text.x = element_text(angle=45, vjust=1, hjust=1)) + 
  labs(size=expression(1/sqrt(SD)),color="Mean Viability\nat Cmax")
ggsave("CmaxViability_heatmap_orderMeanViability_sizeSD-0.5Power.pdf", width=4,height=4)
mean.viability.cmax$size <- mean.viability.cmax$sdViability^-0.5

ggplot(mean.viability.cmax, aes(x=PDX, y=reorder(Drug, -meanViability),
                                color=meanViability, size=sdViability^-0.5)) + 
  geom_point() + 
  facet_wrap(.~timeD)+theme_classic() + scale_size_continuous(breaks=c(0,5,10,20,30,40)^-0.5, 
                                                              labels=c(0,5,10,20,30,40)) + 
  scale_color_gradient(high="grey",low="red",limits=c(0,ceiling(max(mean.viability.cmax$meanViability)))) + 
  theme(axis.title = element_blank(), axis.text.x = element_text(angle=45, vjust=1, hjust=1)) + 
  labs(size="SD",color="Mean Viability\nat Cmax")
ggsave("CmaxViability_heatmap_orderMeanViability_sizeInverseSqrtSD_v2.pdf", width=4,height=4)

ggplot(mean.viability.cmax, aes(x=PDX, y=reorder(Drug, -meanViability),
                                color=meanViability, size=sdViability^-0.5)) + 
  geom_point() + 
  facet_wrap(.~timeD)+theme_classic() + scale_size_continuous(breaks=c(2.5,5,10,20)^-0.5, 
                                                              labels=c(2.5,5,10,20)) + 
  scale_color_gradient(high="grey",low="red",limits=c(0,ceiling(max(mean.viability.cmax$meanViability)))) + 
  theme(axis.title = element_blank(), axis.text.x = element_text(angle=45, vjust=1, hjust=1)) + 
  labs(size="SD",color="Mean Viability\nat Cmax")
ggsave("CmaxViability_heatmap_orderMeanViability_sizeInverseSqrtSD_v3.pdf", width=4,height=4)

mean.viability.cmax[mean.viability.cmax$Drug=="Vorinostat",]$Drug <- "vorinostat"
ggplot(mean.viability.cmax, aes(x=PDX, y=reorder(Drug, -meanViability),
                                color=meanViability, size=sdViability^-0.5)) + 
  geom_point() + 
  facet_wrap(.~timeD)+theme_classic() + scale_size_continuous(breaks=c(2.5,5,10,20)^-0.5, 
                                                              labels=c(2.5,5,10,20)) + 
  scale_color_gradient(high="grey",low="red",limits=c(0,ceiling(max(mean.viability.cmax$meanViability)))) + 
  theme(axis.title = element_blank(), axis.text.x = element_text(angle=45, vjust=1, hjust=1)) + 
  labs(size="SD",color="Mean Viability\nat Cmax")
ggsave("CmaxViability_heatmap_orderMeanViability_sizeInverseSqrtSD_v4.pdf", width=4,height=4)

ggplot(mean.viability.cmax, aes(x=PDX, y=reorder(Drug, -meanViability),
                                color=meanViability, size=meanViability/sdViability)) + 
  geom_point() + facet_wrap(.~timeD)+theme_classic() +  
  scale_color_gradient(high="grey",low="red",limits=c(0,ceiling(max(mean.viability.cmax$meanViability)))) + 
  theme(axis.title = element_blank(), axis.text.x = element_text(angle=45, vjust=1, hjust=1)) + 
  labs(size="Mean/SD",color="Mean Viability\nat Cmax")
ggsave("CmaxViability_heatmap_orderMeanViability_sizeInverseFracSD.pdf", width=4,height=4)

# AUC
auc <- viability[viability$dose_response_metric=="fit_auc",]
quantiles <- quantile(auc$dose_response_value)
auc$Quantile <- NA
for (q in names(quantiles)) {
  if (any(auc$dose_response_value >= quantiles[[q]])){
    auc[auc$dose_response_value >= quantiles[[q]],]$Quantile <- q 
  }
}
auc$Quantile <- factor(auc$Quantile, names(quantiles))
#auc$time <- paste0(auc$time, "h")
#auc$time <- factor(auc$time, levels=c("48h", "120h"))
mean.auc <- plyr::ddply(auc, .(improve_drug_id), summarize,
                        auc=mean(dose_response_value, na.rm=TRUE))

# r2
r2 <- viability[viability$dose_response_metric=="fit_r2",]
quantiles <- quantile(r2$dose_response_value)
r2$Quantile <- NA
for (q in names(quantiles)) {
  if (any(r2$dose_response_value >= quantiles[[q]])){
    r2[r2$dose_response_value >= quantiles[[q]],]$Quantile <- q 
  }
}
r2$Quantile <- factor(r2$Quantile, names(quantiles))

# AUC with r2 fill
auc.r2 <- merge(auc, r2, by=c("improve_drug_id","improve_sample_id", "time"))
drugOrder <- names(drug.info)

ggplot2::ggplot(auc.r2, 
                aes(x=as.factor(time), y=dose_response_value.x, color=dose_response_value.y, 
                    shape=improve_sample_id)) + 
  geom_point() + theme_classic(base_size=12) + #scale_x_discrete() +
  labs(y="Sensitivity Score (AUC)", color=expression(paste("Fit R"^"2")), shape="MPNST") + 
  ggrepel::geom_text_repel(aes(label=improve_drug_id))+
  theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=45, vjust=1, hjust=1))
ggsave(paste0("fit_auc_r2_dotPlot_",Sys.Date(),".pdf"), width=5, height=4)

auc.r2$timeD <- paste0(auc.r2$time/24,"d")
ggplot2::ggplot(auc.r2, 
                aes(x=improve_drug_id, y=dose_response_value.x, fill=dose_response_value.y)) + 
  geom_bar(stat="identity") + theme_classic(base_size=12) + facet_grid(timeD~improve_sample_id) +
  labs(y="Sensitivity Score (AUC)", fill=expression(paste("Fit R"^"2"))) + scale_x_discrete(limits=names(drug.info))+
  theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=45, vjust=1, hjust=1))
ggsave(paste0("fit_auc_r2_barPlot_",Sys.Date(),".pdf"), width=16, height=4) # was width 10 with 3 MPNSTs

# AUC histogram and density plots
hist(r2$dose_response_value, xlab=expression(paste(R^2)))
ggplot(r2, aes(x=dose_response_value, color=improve_sample_id)) + geom_density() + facet_wrap(.~paste0(time/24,"d"))+theme_classic()+xlab(expression(paste(R^2)))
ggsave("R2_density_timeMPNST.pdf",width=8, height=4)

ggplot(r2, aes(x=dose_response_value, fill=improve_sample_id)) + geom_histogram() + facet_wrap(.~paste0(time/24,"d"))+theme_classic()+labs(x=expression(paste(R^2)), fill="MPNST")
ggsave("R2_histogram_timeMPNST.pdf",width=8, height=4)

ggplot(r2, aes(x=dose_response_value)) + geom_histogram() + facet_grid(improve_sample_id~paste0(time/24,"d"))+theme_classic()+xlab(expression(paste(R^2)))
ggsave("R2_histogram_timeMPNSTfacet.pdf",width=8, height=6)

# AUC heatmap
ggplot(auc.r2, aes(x=improve_sample_id, y=reorder(improve_drug_id, -dose_response_value.x), fill=dose_response_value.x)) + geom_tile(stat="identity") + 
  facet_wrap(.~timeD)+theme_classic() + scale_fill_gradient(high="grey",low="red",limits=c(0,1)) + 
  theme(axis.title = element_blank(), axis.text.x = element_text(angle=45, vjust=1, hjust=1)) + labs(fill="AUC")
ggsave("fit_AUC_heatmap.pdf", width=4,height=4)

# r2 > 25% quantile (0.299)
ggplot(auc.r2[auc.r2$Quantile.y != "0%",], aes(x=improve_sample_id, y=reorder(improve_drug_id, -dose_response_value.x), fill=dose_response_value.x)) + geom_tile(stat="identity") + 
  facet_wrap(.~timeD)+theme_classic() + scale_fill_gradient(high="grey",low="red",limits=c(0,1)) + 
  theme(axis.title = element_blank(), axis.text.x = element_text(angle=45, vjust=1, hjust=1)) + labs(fill="AUC")
ggsave("fit_AUC_heatmap_r2PercentileMin25.pdf", width=4,height=4)

# r2 > 50% quantile (0.748)
ggplot(auc.r2[!(auc.r2$Quantile.y %in% c("0%","25%")),], aes(x=improve_sample_id, y=reorder(improve_drug_id, -dose_response_value.x), fill=dose_response_value.x)) + geom_tile(stat="identity") + 
  facet_wrap(.~timeD)+theme_classic() + scale_fill_gradient(high="grey",low="red",limits=c(0,1)) + 
  theme(axis.title = element_blank(), axis.text.x = element_text(angle=45, vjust=1, hjust=1)) + labs(fill="AUC")
ggsave("fit_AUC_heatmap_r2PercentileMin50.pdf", width=4,height=4)

# 75% is 0.898, 100% is 0.983

# #### AUC vs musyc logAlpha ####
# musyc <- read.csv("Deconvolved_musyc_20250804.csv")
# #logAlpha12 is drug1's effect on potency of drug2 so we want to correlate drug1 AUC vs. logAlpha12
# musyc.auc1 <- merge(auc[,c("improve_drug_id","dose_response_value")], musyc[,c("drug1","log_alpha12")], 
#                    by.x="improve_drug_id", by.y="drug1")
# colnames(musyc.auc1) <- c("drug","AUC","log_alpha")
# musyc.auc2 <- merge(auc[,c("improve_drug_id","dose_response_value")], musyc[,c("drug2","log_alpha21")], 
#                     by.x="improve_drug_id", by.y="drug2")
# colnames(musyc.auc2) <- c("drug","AUC","log_alpha")
# musyc.auc <- rbind(musyc.auc1, musyc.auc2) # 3388 across 11 drugs
# cor.test(musyc.auc$AUC, musyc.auc$log_alpha) # weak: r=0.07344482, p=1.876E-5
# length(unique(musyc.auc$drug))
# 
# # try filtering for r2 above median (50+% quartile)
# musyc.auc1 <- merge(auc.r2[auc.r2$Quantile.y %in% c("50%","100%"),
#                            c("improve_drug_id","dose_response_value.x")], musyc[,c("drug1","log_alpha12")], 
#                     by.x="improve_drug_id", by.y="drug1")
# colnames(musyc.auc1) <- c("drug","AUC","log_alpha")
# musyc.auc2 <- merge(auc.r2[auc.r2$Quantile.y %in% c("50%","100%"),
#                            c("improve_drug_id","dose_response_value.x")], musyc[,c("drug2","log_alpha21")], 
#                     by.x="improve_drug_id", by.y="drug2")
# colnames(musyc.auc2) <- c("drug","AUC","log_alpha")
# musyc.auc <- rbind(musyc.auc1, musyc.auc2) # 740 across 9 drugs
# cor.test(musyc.auc$AUC, musyc.auc$log_alpha) # weak: r=0.1193637, p=1.141E-3
# length(unique(musyc.auc$drug))

#### repeat using auc instead of fit_auc ####
## plot
# AUC
auc <- viability[viability$dose_response_metric=="auc",]
quantiles <- quantile(auc$dose_response_value)
auc$Quantile <- NA
for (q in names(quantiles)) {
  if (any(auc$dose_response_value >= quantiles[[q]])){
    auc[auc$dose_response_value >= quantiles[[q]],]$Quantile <- q 
  }
}
auc$Quantile <- factor(auc$Quantile, names(quantiles))
#auc$time <- paste0(auc$time, "h")
#auc$time <- factor(auc$time, levels=c("48h", "120h"))
mean.auc <- plyr::ddply(auc, .(improve_drug_id), summarize,
                        auc=mean(dose_response_value, na.rm=TRUE))

# r2
r2 <- viability[viability$dose_response_metric=="fit_r2",]
quantiles <- quantile(r2$dose_response_value)
r2$Quantile <- NA
for (q in names(quantiles)) {
  if (any(r2$dose_response_value >= quantiles[[q]])){
    r2[r2$dose_response_value >= quantiles[[q]],]$Quantile <- q 
  }
}
r2$Quantile <- factor(r2$Quantile, names(quantiles))

# AUC with r2 fill
auc.r2 <- merge(auc, r2, by=c("improve_drug_id","improve_sample_id", "time"))
drugOrder <- names(drug.info)

ggplot2::ggplot(auc.r2, 
                aes(x=as.factor(time), y=dose_response_value.x, color=dose_response_value.y, 
                    shape=improve_sample_id)) + 
  geom_point() + theme_classic(base_size=12) + #scale_x_discrete() +
  labs(y="Sensitivity Score (AUC)", color=expression(paste("Fit R"^"2")), shape="MPNST") + 
  ggrepel::geom_text_repel(aes(label=improve_drug_id))+
  theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=45, vjust=1, hjust=1))
ggsave(paste0("auc_r2_dotPlot_",Sys.Date(),".pdf"), width=5, height=4)

auc.r2$timeD <- paste0(auc.r2$time/24,"d")
ggplot2::ggplot(auc.r2, 
                aes(x=improve_drug_id, y=dose_response_value.x, fill=dose_response_value.y)) + 
  geom_bar(stat="identity") + theme_classic(base_size=12) + facet_grid(timeD~improve_sample_id) +
  labs(y="Sensitivity Score (AUC)", fill=expression(paste("Fit R"^"2"))) + scale_x_discrete(limits=names(drug.info))+
  theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=45, vjust=1, hjust=1))
ggsave(paste0("auc_r2_barPlot_",Sys.Date(),".pdf"), width=16, height=4) # was width 10 with 3 MPNSTs

# AUC heatmap
ggplot(auc.r2, aes(x=improve_sample_id, y=reorder(improve_drug_id, -dose_response_value.x), fill=dose_response_value.x)) + geom_tile(stat="identity") + 
  facet_wrap(.~timeD)+theme_classic() + scale_fill_gradient(high="grey",low="red",limits=c(0,1)) + 
  theme(axis.title = element_blank(), axis.text.x = element_text(angle=45, vjust=1, hjust=1)) + labs(fill="AUC")
ggsave("AUC_heatmap.pdf", width=4,height=4)

# r2 > 25% quantile (0.530800)
ggplot(auc.r2[auc.r2$Quantile.y != "0%",], aes(x=improve_sample_id, y=reorder(improve_drug_id, -dose_response_value.x), fill=dose_response_value.x)) + geom_tile(stat="identity") + 
  facet_wrap(.~timeD)+theme_classic() + scale_fill_gradient(high="grey",low="red",limits=c(0,1)) + 
  theme(axis.title = element_blank(), axis.text.x = element_text(angle=45, vjust=1, hjust=1)) + labs(fill="AUC")
ggsave("AUC_heatmap_r2PercentileMin25.pdf", width=4,height=4)

# r2 > 50% quantile (0.799950)
ggplot(auc.r2[!(auc.r2$Quantile.y %in% c("0%","25%")),], aes(x=improve_sample_id, y=reorder(improve_drug_id, -dose_response_value.x), fill=dose_response_value.x)) + geom_tile(stat="identity") + 
  facet_wrap(.~timeD)+theme_classic() + scale_fill_gradient(high="grey",low="red",limits=c(0,1)) + 
  theme(axis.title = element_blank(), axis.text.x = element_text(angle=45, vjust=1, hjust=1)) + labs(fill="AUC")
ggsave("AUC_heatmap_r2PercentileMin50.pdf", width=4,height=4)

filtered.auc <- auc.r2[!(auc.r2$Quantile.y %in% c("0%","25%")),]
mean.filt.auc <- plyr::ddply(filtered.auc, .(improve_drug_id), summarize,
                             meanAUC = mean(dose_response_value.x),
                             medianAUC = mean(dose_response_value.x),
                             maxAUC = max(dose_response_value.x),
                             minAUC = min(dose_response_value.x))
ggplot(auc.r2, aes(x=improve_sample_id, y=improve_drug_id,
                   color=dose_response_value.x, size=dose_response_value.y)) + 
  geom_point() + scale_y_discrete(limits=mean.filt.auc[order(-mean.filt.auc$meanAUC),]$improve_drug_id)+
  facet_wrap(.~timeD)+theme_classic() + scale_color_gradient(high="grey",low="red",limits=c(0,1)) + 
  theme(axis.title = element_blank(), axis.text.x = element_text(angle=45, vjust=1, hjust=1)) + labs(size=expression(paste(R^2)),color="AUC")
ggsave("AUC_heatmap_orderMeanR2PercentileMin50_sizeR2.pdf", width=4,height=4)

ggplot(auc.r2, aes(x=improve_sample_id, y=reorder(improve_drug_id, -dose_response_value.x),
                   color=dose_response_value.x, size=dose_response_value.y)) + 
  geom_point() + 
  facet_wrap(.~timeD)+theme_classic() + scale_color_gradient(high="grey",low="red",limits=c(0,1)) + 
  theme(axis.title = element_blank(), axis.text.x = element_text(angle=45, vjust=1, hjust=1)) + labs(size=expression(paste(R^2)),color="AUC")
ggsave("AUC_heatmap_orderMeanAUC_sizeR2.pdf", width=4,height=4)

# 75% is 0.922925, 100% is 0.996100
# 
# #### AUC vs musyc logAlpha ####
# musyc <- read.csv("Deconvolved_musyc_20250804.csv")
# #logAlpha12 is drug1's effect on potency of drug2 so we want to correlate drug1 AUC vs. logAlpha12
# musyc.auc1 <- merge(auc[,c("improve_drug_id","dose_response_value")], musyc[,c("drug1","log_alpha12")], 
#                     by.x="improve_drug_id", by.y="drug1")
# colnames(musyc.auc1) <- c("drug","AUC","log_alpha")
# musyc.auc2 <- merge(auc[,c("improve_drug_id","dose_response_value")], musyc[,c("drug2","log_alpha21")], 
#                     by.x="improve_drug_id", by.y="drug2")
# colnames(musyc.auc2) <- c("drug","AUC","log_alpha")
# musyc.auc <- rbind(musyc.auc1, musyc.auc2) # 3388 across 11 drugs
# cor.test(musyc.auc$AUC, musyc.auc$log_alpha) # no: r = -0.02678685, p = 0.119
# length(unique(musyc.auc$drug))
# 
# # try filtering for r2 above 25% quartile
# musyc.auc1 <- merge(auc.r2[auc.r2$Quantile.y %in% c("25%","50%","100%"),
#                            c("improve_drug_id","dose_response_value.x")], musyc[,c("drug1","log_alpha12")], 
#                     by.x="improve_drug_id", by.y="drug1")
# colnames(musyc.auc1) <- c("drug","AUC","log_alpha")
# musyc.auc2 <- merge(auc.r2[auc.r2$Quantile.y %in% c("25%","50%","100%"),
#                            c("improve_drug_id","dose_response_value.x")], musyc[,c("drug2","log_alpha21")], 
#                     by.x="improve_drug_id", by.y="drug2")
# colnames(musyc.auc2) <- c("drug","AUC","log_alpha")
# musyc.auc <- rbind(musyc.auc1, musyc.auc2) # 1760 across 11 drugs
# cor.test(musyc.auc$AUC, musyc.auc$log_alpha) # weak: r = -0.0.05054397, p = 3.398E-2
# length(unique(musyc.auc$drug))
# 
# # try filtering for r2 above median (50+% quartile)
# musyc.auc1 <- merge(auc.r2[auc.r2$Quantile.y %in% c("50%","100%"),
#                            c("improve_drug_id","dose_response_value.x")], musyc[,c("drug1","log_alpha12")], 
#                     by.x="improve_drug_id", by.y="drug1")
# colnames(musyc.auc1) <- c("drug","AUC","log_alpha")
# musyc.auc2 <- merge(auc.r2[auc.r2$Quantile.y %in% c("50%","100%"),
#                            c("improve_drug_id","dose_response_value.x")], musyc[,c("drug2","log_alpha21")], 
#                     by.x="improve_drug_id", by.y="drug2")
# colnames(musyc.auc2) <- c("drug","AUC","log_alpha")
# musyc.auc <- rbind(musyc.auc1, musyc.auc2) # 740 across 9 drugs
# cor.test(musyc.auc$AUC, musyc.auc$log_alpha) # weak: r = -0.1217981, p = 9.001E-4
# length(unique(musyc.auc$drug))
