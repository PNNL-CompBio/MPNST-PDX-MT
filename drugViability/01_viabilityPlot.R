# viability plot
library(plyr);library(dplyr);library(synapser);library(ggplot2)
synapser::synLogin()
base.path <- "~/OneDrive - PNNL/Documents/GitHub/MPNST-PDX-MT/drugViability"
setwd(base.path)

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

## plot
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
ggsave("fit_auc_r2_dotPlot_2025-05-21.pdf", width=5, height=4)

auc.r2$timeD <- paste0(auc.r2$time/24,"d")
ggplot2::ggplot(auc.r2, 
                aes(x=improve_drug_id, y=dose_response_value.x, fill=dose_response_value.y)) + 
  geom_bar(stat="identity") + theme_classic(base_size=12) + facet_grid(timeD~improve_sample_id) +
  labs(y="Sensitivity Score (AUC)", fill=expression(paste("Fit R"^"2"))) + scale_x_discrete(limits=names(drug.info))+
  theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=45, vjust=1, hjust=1))
ggsave("fit_auc_r2_barPlot_2025-05-21.pdf", width=10, height=4)
