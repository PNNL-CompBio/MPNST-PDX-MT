# viability plot
library(plyr);library(dplyr);library(synapser);library(ggplot2)
synapser::synLogin()
#base.path <- "~/OneDrive - PNNL/Documents/GitHub/MPNST-PDX-MT/EXACT-pipeline/drugViability"
#setwd(base.path)
setwd("./")
# run 00_createCurveStats.py
library(reticulate)

#source_python("00_createCurveStats.py")

drug.info <- list("palbociclib" = "CDK inhibitor", "ribociclib" = "CDK inhibitor",
                  "trabectedin" = "Chemotherapy", "ifosfamide" = "DNA alkylating agent",
                  "decitabine" = "DNMT inhibitor", "vorinostat" = "HDAC inhibitor",
                  "mirdametinib" = "MEK inhibitor", "selumetinib" = "MEK inhibitor",
                  "trametinib" = "MEK inhibitor", "capmatinib" = "MET inhibitor",
                  "olaparib" = "PARP inhibitor", "RMC4630" = "SHP2 inhibitor",
                  "TNO155" = "SHP2 inhibitor", "doxorubicin" = "TOP inhibitor",
                  "irinotecan" = "TOP inhibitor", "verteporfin" = "YAP inhibitor")

## load data
viability <- read.table(synapser::synGet("syn65941820")$path, sep = "\t", header = TRUE)
viability <- viability[viability$improve_drug_id != "irinotecan",]
viability <- viability[viability$improve_drug_id != "avutometinib",]
viability <- viability[viability$improve_drug_id != "defactinib",]
viability <- viability[viability$improve_drug_id != "nazartinib",]
viability[viability$improve_drug_id == "SN38",]$improve_drug_id <- "SN-38"
viability.cmax <- readxl::read_excel(synapser::synGet("syn70107593")$path)


## plot
# viability at Cmax
mean.viability.cmax <- plyr::ddply(viability.cmax, .(Drug, PDX, Hour), summarize,
                                   meanViability = mean(`Viability at Cmax`, na.rm=TRUE),
                                   sdViability = sd(`Viability at Cmax`, na.rm=TRUE))
mean.viability.cmax$timeD <- paste0(as.numeric(mean.viability.cmax$Hour)/24,"d")


if ('Vorinostat'%in%mean.viability.cmax$Drug)
  mean.viability.cmax[mean.viability.cmax$Drug=="Vorinostat",]$Drug <- "vorinostat"
mean.viability.cmax[mean.viability.cmax$Drug=="SN38",]$Drug <- "SN-38"
if ('JH2002' %in% mean.viability.cmax$PDX)
  mean.viability.cmax[mean.viability.cmax$PDX=="JH2002",]$PDX <- "JH-2-002"
if ('JH2079c' %in% mean.viability.cmax$PDX)
  mean.viability.cmax[mean.viability.cmax$PDX=="JH2079c",]$PDX <- "JH-2-079c"
if('MN2' %in% mean.viability.cmax$PDX)
  mean.viability.cmax[mean.viability.cmax$PDX=="MN2",]$PDX <- "MN-2"
if('MN4' %in% mean.viability.cmax$PDX)
  mean.viability.cmax[mean.viability.cmax$PDX=="MN4",]$PDX <- "MN-4"
if('WU225' %in% mean.viability.cmax$PDX)
  mean.viability.cmax[mean.viability.cmax$PDX=="WU225",]$PDX <- "WU-225"
if('WU356' %in% mean.viability.cmax$PDX)
  mean.viability.cmax[mean.viability.cmax$PDX=="WU356",]$PDX <- "WU-356"

ggplot(mean.viability.cmax, aes(x=PDX, y=reorder(Drug, -meanViability),
                                color=meanViability, size=sdViability^-0.5)) +
  geom_point() +
  facet_wrap(.~timeD)+theme_classic() + scale_size_continuous(breaks=c(2.5,5,10,20)^-0.5,
                                                              labels=c(2.5,5,10,20)) +
  scale_color_gradient(high="grey",low="red",limits=c(0,ceiling(max(mean.viability.cmax$meanViability)))) +
  theme(axis.title = element_blank(), axis.text.x = element_text(angle=45, vjust=1, hjust=1)) +
  labs(size="SD",color="Mean\nViability\nat Cmax")
ggsave("fig_2B_CmaxViability_heatmap_orderMeanViability_sizeInverseSqrtSD_final.pdf", width=4,height=4)

ggplot(mean.viability.cmax, aes(x=PDX, y=reorder(Drug, -meanViability, FUN=median),
                                color=meanViability, size=sdViability^-0.5)) +
  geom_point() +
  facet_wrap(.~timeD)+theme_classic() + scale_size_continuous(breaks=c(2.5,5,10,20)^-0.5,
                                                              labels=c(2.5,5,10,20)) +
  scale_color_gradient(high="grey",low="red",limits=c(0,ceiling(max(mean.viability.cmax$meanViability)))) +
  theme(axis.title = element_blank(), axis.text.x = element_text(angle = 45, vjust=1, hjust=1)) +
  labs(size="SD",color="Mean\nViability\nat Cmax")
#ggsave("CmaxViability_heatmap_orderMedianViability_sizeInverseSqrtSD_final.pdf", width=4,height=4)

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
auc.r2 <- rbind(auc,r2) |>#merge(auc, r2, by=c("improve_drug_id","improve_sample_id", "time"))
    dplyr::select(-Quantile)|>
  tidyr::pivot_wider(names_from='dose_response_metric',values_from='dose_response_value')
auc.r2$timeD = paste0(as.numeric(auc.r2$time)/24,"d")
ggplot(auc.r2, aes(x=improve_sample_id, y=reorder(improve_drug_id, -auc),
                   color=auc, size=fit_r2)) +
  geom_point() +
  facet_wrap(.~timeD)+theme_classic() + scale_color_gradient(high="grey",low="red",limits=c(0,1)) +
  theme(axis.title = element_blank(), axis.text.x = element_text(angle=45, vjust=1, hjust=1)) + labs(size=expression(paste(R^2)),color="AUC")
#ggsave(paste0("AUC_heatmap_orderMeanAUC_sizeR2_",Sys.Date(),".pdf"), width=4,height=4)

ggplot(auc.r2, aes(x=improve_sample_id, y=reorder(improve_drug_id, -auc, FUN=median),
                   color=auc, size=fit_r2)) +
  geom_point() +
  facet_wrap(.~timeD)+theme_classic() + scale_color_gradient(high="grey",low="red",limits=c(0,1)) +
  theme(axis.title = element_blank(), axis.text.x = element_text(angle=45, vjust=1, hjust=1)) + labs(size=expression(paste(R^2)),color="AUC")
ggsave(paste0("fig_2A_AUC_heatmap_orderMedianAUC_sizeR2_",Sys.Date(),".pdf"), width=4,height=4)
