  ##load drug data for DMEA
#moving to separate file for use by other tools


##function to lod mRNA
library(synapser)
# load helper functions
get_CCLE_RNA <- function(hgnc = FALSE) {
  message("Loading adherent CCLE RNA-seq data version 19Q4")
  if (hgnc) {


    RNA.first200 <- readRDS(synGet('syn74935645')$path)##"Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_approved_symbols_only_1-200.Rbin")
    RNA.rest <- readRDS(synGet('syn74935650')$path)#"Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_approved_symbols_only_201-327.Rbin")
    RNA.df <- rbind(RNA.first200, RNA.rest)
  } else {


    load(synGet('syn74954979')$path)#"Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_1-200.Rbin")
    load(synGet('syn74954995')$path)#"Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_201-327.Rbin")
    RNA.df <- rbind(RNA.first200, RNA.rest)
  }
  return(RNA.df)
}



# PRISM drug sensitivity (i.e., AUC)
prism <- read.csv(synGet('syn74955119')$path)
prism$X <- NULL

# PRISM drug mechanism of action annotations
gmt.drug <- GSA::GSA.read.gmt(synGet('syn74955107')$path)
