# create supp table for diffexp
library(synapser);library(openxlsx)
synapser::synLogin()
# read in diffexp
individual <- read.csv(synapser::synGet("syn69910880")$path)
overall <- read.csv(synapser::synGet("syn64397743")$path)

readme <- data.frame("Sheet" = c("individualSamples","overall"), 
                     "Description" = c("Differential expression for each MPNST at each time for each drug relative to DMSO controls.",
                                       "Differential expression at each time for each drug relative to DMSO controls."))

de <- list("README" = readme, "individualSamples" = individual, "overall" = overall)
temp.path <- "/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/MPNST/Synergy"
openxlsx::write.xlsx(de, 
                     file=file.path(temp.path,
                                    paste0("SupplementaryTable",2,"_","differentialExpression",".xlsx")), 
                     rowNames=FALSE)
