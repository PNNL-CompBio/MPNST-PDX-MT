# Project
MPNST-PDX-MT-RNAseq_data
## Batches
1. batch1: Mirda on JH-2-002
2. batch2: All drugs on MN2
## Analysis scripts which read to and from Synapse
### Differential expression analysis
rnaseqSummary.Rmd
#### Pulls: sample counts
#### Uploads: counts crosstab; differential expression

### Transcription factor enrichment analysis
tfEnrichment.Rmd
#### Pulls: counts crosstab; differential expression
#### Uploads: transcription factor enrichment

### Drug mechanism enrichment analysis (DMEA)
1. DMEA.Rmd
#### Pulls: differential expression; public adherent CCLE RNA-seq and PRISM drug sensitivity
#### Uploads: DMEA results based on sensitivity

2. DMEA_L1000.Rmd
#### Pulls: differential expression; public Connectivity Map L1000
#### Uploads: DMEA results based on similarity (i.e., correlations)

3. DMEA_Cmap_L1000.Rmd
#### Pulls: differential expression; public Connectivity Map L1000
#### Uploads: DMEA results based on similarity (i.e., clue.io Cmap L1000 Query)