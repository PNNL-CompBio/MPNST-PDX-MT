# MEK+HDACi differential expression
library(plyr);library(dplyr);library(synapser);
library(ggplot2);library(ComplexHeatmap);library(RColorBrewer); library(ggvenn)
#source("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/helperScripts/compile_mDEG.R")
library(panSEA)

setwd("~/OneDrive - PNNL/Documents/MPNST")
dir.create("MEKi")
setwd("MEKi")
dir.create("Treated")
setwd("Treated")
dir.create("GSE183307_and_GSE262030")
setwd("GSE183307_and_GSE262030")
dir.create("24h")
setwd("24h")
dir.create("padj0.05")
setwd("padj0.05")
dir.create("min2PerGroup")
setwd("min2PerGroup")

#### GSE183307 ####
# load counts table from GEO
url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE183307&format=file&file=GSE183307%5FPartek%5Frawreadscount%5FMEKi%2Ecsv%2Egz"
#path <- paste(urld, "acc=GSE183307", "file=GSE183307_Partek_rawreadscount_MEKi.csv.gz", sep="&");
tbl <- as.matrix(data.table::fread(url, header=T, colClasses="integer"), rownames="Gene Symbol") # 237372

# create metadata
Sample <- colnames(tbl)
factor.info <- data.frame(Sample, Treatment="DMSO", Time="8h")
factor.info[grepl("24h",factor.info$Sample),]$Time <- "24h"
factor.info[grepl("MEKi",factor.info$Sample),]$Treatment <- "MEKi"
rownames(factor.info) <- factor.info$Sample
factor.info$Sample <- NULL
factor.info$Treatment <- factor(factor.info$Treatment, levels=c("MEKi","DMSO"))
factor.info$Time <- factor(factor.info$Time, levels=c("8h","24h"))
factor.info <- factor.info[factor.info$Time == "24h",]
factor.info$Time <- NULL
tbl <- tbl[,rownames(factor.info)]

# pre-filter low count genes
# keep genes with at least N+1 counts > 10, where N = size of smallest group (3)
keep <- rowSums( tbl >= 10 ) > 3
tbl <- tbl[keep, ]

# # only keep genes quantified in at least half of samples
# missing <- rowSums(tbl == 0)
# tbl <- tbl[missing < ncol(tbl)/2, ]

# further, require genes to be quantified in at least 2 samples per group
missingDMSO <- rowSums(tbl[grepl("DMSO", colnames(tbl)),] == 0)
tbl <- tbl[missingDMSO < 2, ]
missingDrug <- rowSums(tbl[grepl("MEKi", colnames(tbl)),] == 0)
tbl <- tbl[missingDrug < 2, ]

# notice COL1A1 is very abundant

# run diffexp
tbl <- tbl[!duplicated(row.names(tbl)),]
eset <- Biobase::ExpressionSet(tbl)
ex <- Biobase::exprs(eset)
qx <- as.numeric(quantile(ex, c(0, 0.25, 0.5, 
                                0.75, 0.99, 1), na.rm = TRUE))
LogC <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && 
                            qx[2] > 0)
if (LogC) {
  ex[which(ex <= 0)] <- NaN
  Biobase::exprs(eset) <- log2(ex)
}
eset$group <- factor.info[, 1]
if (ncol(factor.info) > 1) {
  eset$batch <- factor.info[, 2]
  levels(eset$batch) <- levels(factor.info[, 
                                           2])
  design <- stats::model.matrix(~group + batch + 
                                  0, eset)
  colnames(design)[1:2] <- substr(colnames(design)[1:2], 
                                  6, nchar(colnames(design)[1:2]))
} else {
  design <- stats::model.matrix(~group + 0, 
                                eset)
  colnames(design) <- substr(colnames(design), 
                             6, nchar(colnames(design)))
}
fit <- limma::lmFit(eset, design)
cts <- paste(levels(factor.info[, 1]), collapse = "-")
cont.matrix <- limma::makeContrasts(contrasts = cts, 
                                    levels = design)
fit <- limma::contrasts.fit(fit, cont.matrix)
fit <- limma::eBayes(fit, 0.01)
deg <- limma::topTable(fit, adjust = "fdr", 
                                   number = nrow(fit))
deg$Gene <- rownames(deg)
colnames(deg)[1] <- "Log2FC"
deg <- deg[, c("Gene", colnames(deg)[1:(ncol(deg) - 1)])]
deg$Log2Transformed <- LogC
write.csv(deg, "GSE183307_diffexp.csv", row.names = FALSE)

#### GSE262030 ####
# load counts table from GEO
urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
path <- paste(urld, "acc=GSE262030", "file=GSE262030_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep="&");
tbl <- as.matrix(data.table::fread(path, header=T, colClasses="integer"), rownames="GeneID")

# load gene annotations 
apath <- paste(urld, "type=rnaseq_counts", "file=Human.GRCh38.p13.annot.tsv.gz", sep="&")
annot <- data.table::fread(apath, header=T, quote="", stringsAsFactors=F, data.table=F)
rownames(annot) <- annot$GeneID

# switch to gene symbols
annot <- annot[!is.na(annot$Symbol) & annot$GeneID %in% rownames(tbl),]
tbl <- tbl[rownames(tbl) %in% annot$GeneID,]
rownames(tbl) <- annot[rownames(tbl),]$Symbol

# sample selection
gsms <- "XXXXXXXXXXX11XXXXXXXXXXXXXXXXXXXXXXXXXXX000XXX"
sml <- strsplit(gsms, split="")[[1]]

# filter out excluded samples (marked as "X")
sel <- which(sml != "X")
sml <- sml[sel]
tbl <- tbl[ ,sel]

# group membership for samples
gs <- factor(sml)
groups <- make.names(c("Selumetinib","DMSO"))
levels(gs) <- groups
sample_info <- data.frame(Group = gs, row.names = colnames(tbl))

# pre-filter low count genes
# keep genes with at least N+1 counts > 10, where N = size of smallest group
keep <- rowSums( tbl >= 10 ) > min(table(gs))
tbl <- tbl[keep, ]

# create metadata
Sample <- colnames(tbl)
Treatment <- c("DMSO","DMSO","Selumetinib","Selumetinib","Selumetinib")
factor.info <- data.frame(Sample, Treatment)
rownames(factor.info) <- factor.info$Sample
factor.info$Sample <- NULL
factor.info$Treatment <- factor(factor.info$Treatment, levels=c("Selumetinib","DMSO"))

# further, require genes to be quantified in at least 2 samples per group
DMSOSamples <- rownames(factor.info)[factor.info$Treatment == "DMSO"]
missingDMSO <- rowSums(tbl[,DMSOSamples] == 0)
tbl <- tbl[missingDMSO < 2, ]
DrugSamples <- rownames(factor.info)[factor.info$Treatment == "Selumetinib"]
missingDrug <- rowSums(tbl[,DrugSamples] == 0)
tbl <- tbl[missingDrug < 2, ]

# run diffexp
tbl <- tbl[!duplicated(row.names(tbl)),]
eset <- Biobase::ExpressionSet(tbl)
ex <- Biobase::exprs(eset)
qx <- as.numeric(quantile(ex, c(0, 0.25, 0.5, 
                                0.75, 0.99, 1), na.rm = TRUE))
LogC <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && 
                            qx[2] > 0)
if (LogC) {
  ex[which(ex <= 0)] <- NaN
  Biobase::exprs(eset) <- log2(ex)
}
eset$group <- factor.info[, 1]
design <- stats::model.matrix(~group +  
                                0, eset)
colnames(design)[1:2] <- substr(colnames(design)[1:2], 
                                6, nchar(colnames(design)[1:2]))
fit <- limma::lmFit(eset, design)
cts <- paste(levels(factor.info[, 1]), collapse = "-")
cont.matrix <- limma::makeContrasts(contrasts = cts, 
                                    levels = design)
fit <- limma::contrasts.fit(fit, cont.matrix)
fit <- limma::eBayes(fit, 0.01)
deg <- limma::topTable(fit, adjust = "fdr", 
                       number = nrow(fit))
deg$Gene <- rownames(deg)
colnames(deg)[1] <- "Log2FC"
deg <- deg[, c("Gene", colnames(deg)[1:(ncol(deg) - 1)])]
deg$Log2Transformed <- LogC
write.csv(deg, "GSE262030_diffexp.csv", row.names = FALSE)
