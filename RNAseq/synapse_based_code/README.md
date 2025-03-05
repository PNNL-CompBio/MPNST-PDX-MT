## RNASeq Analysis pipeline

This folder contains a list of steps to analyze the RNASeq from Synapse

### 00_rnaseqSummary.Rmd
Collects basic data from synapse and collapses into PCA, then does diferential expression

### 01_tfEnrichment.Rmd
Runs and analyzes decoupleR tf enrichment. 

### 02_priorTFEnrichment.Rmd
Runs and analyzes priori tf enrichment. 

### 03_compareRNASeqAnalyses.Rmd

Here we compare the two TF enrichments to see which is more reliable in summarizing the data and experiment with ways of visualizing.
