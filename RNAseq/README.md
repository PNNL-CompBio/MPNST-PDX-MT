# Project
MPNST-PDX-MT-RNAseq_data
# Batches
1. batch1: Mirda on  MN2
2. batch2: Other drugs on MN2

# Running code (download rmd files, Rfunction.R and config.yml in your local computer; set up config.yml before running)
1. Running Counts_TPM_batch1_2.Rmd 
   * This code is used for checking data "Batch1 and Batch2" dowloaded from synapse. The output pdf is in the synapse syn64612888
   
## Analysis scripts which read to and from Synapse
### Differential expression analysis
rnaseqSummary.Rmd
- Pulls: sample counts (syn52369043)
- Uploads: counts crosstabs (syn64398049 and syn64398051); differential expression (syn64397743)

### Transcription factor enrichment analysis
tfEnrichment.Rmd
- Pulls: counts crosstabs (syn64398049 and syn64398051); differential expression (syn64397743)
- Uploads: transcription factor enrichment (syn64693461)

### Drug mechanism enrichment analysis (DMEA)
1. DMEA.Rmd
- Pulls: differential expression (syn64397743); public adherent CCLE RNA-seq and PRISM drug sensitivity
- Uploads: DMEA results based on sensitivity (syn64618500)

2. DMEA_L1000.Rmd
- Pulls: differential expression (syn64397743); public Connectivity Map L1000
- Uploads: DMEA results based on similarity (i.e., correlations) (syn64660160)

3. DMEA_Cmap_L1000.Rmd
- Pulls: differential expression (syn64397743); public Connectivity Map L1000
- Uploads: DMEA results based on similarity (i.e., clue.io Cmap L1000 Query) (syn64693139)