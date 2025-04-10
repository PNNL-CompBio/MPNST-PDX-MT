# viability plot
library(plyr);library(dplyr);library(synapser);library(ggplot2)
synapser::synLogin()
base.path <- "~/OneDrive - PNNL/Documents/GitHub/MPNST-PDX-MT/drugViability"
setwd(base.path)

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
auc$time <- paste0(auc$time, "h")
auc$time <- factor(auc$time, levels=c("48h", "120h"))
mean.auc <- plyr::ddply(auc, .(improve_drug_id), summarize,
                        auc=mean(dose_response_value, na.rm=TRUE))
drugOrder <- unique(mean.auc[order(mean.auc$auc),]$improve_drug_id)
ggplot2::ggplot(auc, aes(x=improve_drug_id, y=dose_response_value, fill=Quantile)) + 
  geom_bar(stat="identity") + theme_classic(base_size=12) + facet_grid(time ~ improve_sample_id) +
  labs(y="Sensitivity Score (AUC)") + scale_x_discrete(limits=drugOrder)+
  theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=45, vjust=1, hjust=1))
ggsave("fit_AUC_barPlot.pdf", width=9, height=4)

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
r2$time <- paste0(r2$time, "h")
r2$time <- factor(r2$time, levels=c("48h", "120h"))
mean.r2 <- plyr::ddply(r2, .(improve_drug_id), summarize,
                        r2=mean(dose_response_value, na.rm=TRUE))
drugOrder <- unique(mean.r2[order(mean.r2$r2),]$improve_drug_id)
ggplot2::ggplot(r2, aes(x=improve_drug_id, y=dose_response_value, fill=Quantile)) + 
  geom_bar(stat="identity") + theme_classic(base_size=12) + facet_grid(time ~ improve_sample_id) +
  labs(y="Fit Score (r2)") + scale_x_discrete(limits=drugOrder)+
  theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=45, vjust=1, hjust=1))
ggsave("fit_r2_barPlot.pdf", width=9, height=4)

# AUC with r2 fill
auc.r2 <- merge(auc, r2, by=c("improve_drug_id","improve_sample_id", "time"))
drugOrder <- unique(mean.auc[order(mean.auc$auc),]$improve_drug_id)
ggplot2::ggplot(auc.r2, aes(x=improve_drug_id, y=dose_response_value.x, fill=Quantile.y)) + 
  geom_bar(stat="identity") + theme_classic(base_size=12) + facet_grid(time ~ improve_sample_id) +
  labs(y="Sensitivity Score (AUC)", fill="Fit Quantile") + scale_x_discrete(limits=drugOrder)+
  theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=45, vjust=1, hjust=1))
ggsave("fit_auc_r2_barPlot.pdf", width=10, height=4)
