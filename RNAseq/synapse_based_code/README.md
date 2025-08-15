# Project
MPNST-PDX-MT-RNAseq_data
# Batches
1. batch1: Mirda on  MN2
2. batch2: Other drugs on MN2
3. batch3: Other drugs on JH-2-002

# Running code (download rmd files, Rfunction.R and config.yml in your local computer; set up config.yml before running)
1. Running Counts_TPM_batch1_2.Rmd 
   * This code is used for checking data "Batch1 and Batch2" dowloaded from synapse. The output pdf is in the synapse syn64612888

   
## Analysis scripts which read to and from Synapse
### Differential expression analysis
00_rnaseqSummary.Rmd
- Pulls: sample counts (syn52369043)
- Uploads: counts crosstabs (syn64398049 and syn64398051); differential expression (syn64397743)

### Transcription factor activity
01_tfEnrichment.Rmd
- Pulls: counts crosstabs (syn64398049 and syn64398051); differential expression (syn64397743)
- Uploads: transcription factor activity (syn64693461)

### Comparison of differential expression across drugs
02_compareRNASeqAnalyses.Rmd
- Pulls: differential expression (syn64397743)

### Gene set enrichment analysis (GSEA)
03_GSEA.Rmd
- Pulls: differential expression (syn64397743)
- Uploads: enrichment of MSigDB hallmark gene sets (syn68702141)

### Drug mechanism enrichment analysis (DMEA)
04_DMEA.Rmd
- Pulls: differential expression (syn64397743); public adherent CCLE RNA-seq and PRISM drug sensitivity
- Uploads: DMEA results based on sensitivity (syn65672258)

### Evidence for MEK inhibitor combos
05_MEKi_combos.Rmd
- Pulls: PDX mirdametinib + vorinostat tumor size (syn68900596); differential expression (syn64397743); transcription factor activity (syn64693461); GSEA (syn68702141); DMEA (syn65672258)

