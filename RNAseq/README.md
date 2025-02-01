# Project
MPNST-PDX-MT-RNAseq_data
## Batches
1. batch1: Mirda on JH-2-002
2. batch2: All drugs on MN2
## Analysis scripts which read to and from Synapse
### Differential expression analysis
rnaseqSummary.Rmd: pulls counts and uploads differential expression

### Transcription factor enrichment analysis
tfEnrichment.Rmd: pulls differential expression and uploads transcription factor enrichment

### Drug mechanism enrichment analysis
1. DMEA.Rmd: pulls differential expression and public adherent CCLE RNA-seq and PRISM drug sensitivity data and uploads drug rankings based on sensitivity
2. DMEA_L1000.Rmd: pulls differential expression and public Connectivity Map L1000 data to run correlations and uploads drug rankings based on similarity
3. DMEA_Cmap_L1000.Rmd: pulls differential expression and public Connectivity Map L1000 data to run clue.io Query and uploads drug rankings based on similarity