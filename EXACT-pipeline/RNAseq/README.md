# RNASeq EXACT analysis

We use the RNA-sequencing data to identify drugs that target similar pathways. Given the 
dependencies across the analyses we recommend running them in order. 



## Differential expression

For each drug and sample, we calculate the difference between treated and 
DMSO samples. 

### Assemble batches

First we assemble the batches, download into a single file, and handle renaming
and batch correction.

1. batch1: Mirda on  MN2
2. batch2: Other drugs on MN2
3. batch3: Other drugs on JH-2-002
4. batch4: WU-225

All data is downloaded via `downloadCorrectRNASeq.R` and processed through
two separate pipelines:

### Individual differential expression
We  calculate the differential expression across individual 3D-MEDS samples
per patient (MN-2, JH-2-002, and WU-225). Run this before the next step. 

This produces table S2: 00_rnaseqSummary_individual.Rmd
- Pulls: metadata counts (syn64608368,syn64608372,syn66050299 and syn72518652)
- Uploads: counts batch-corrected crosstabs (syn64398049 and syn64398051); 
differential expression for MPNST, time, and drug pairs (syn69910880)

### Combined differential expression
We then combine the measurements across individuals (MN-2, JH-2-002, and WU-225)
to assess differential expression consistent across all samples

This produces another table S2: 00_rnaseqSummary.Rmd
- Pulls: batch-corrected cross tab (syn64398049)
- Uploads:differential expression for time and drug pairs (syn64397743)


### Figure 3A: Comparison of differential expression across drugs
This file counts how many genes are differentially expressed across
time and drugs.

02_compareRNASeqAnalyses.Rmd
- Pulls: differential expression for time and drug pairs (syn64397743)

## TF Target analysis
We then compute the TF active using decoupler

### Figure 3B Transcription factor activity
01_tfEnrichment.Rmd
- Pulls: counts crosstabs (syn64398049 and syn64398051); 
differential expression for time and drug pairs (syn64397743)
- Uploads: transcription factor activity for time and drug pairs (syn69911032)



## Gene set enrichment analysis

We run the GSEA analysis using cancer hallmarks across both individual and 
pooled samples using the `03a_runGSEA.R` script. This will store
the results in individual directories for ploting. 
### Figure 3C: Gene set enrichment analysis (GSEA) across samples

Figure 3B, Table S3: 03_GSEA_Hallmark.Rmd
- Pulls: differential expression for time and drug pairs (syn64397743)
- Uploads: enrichment of MSigDB hallmark gene sets for 
time and drug pairs (syn68702141)


## NOT RUN: 
01_tfEnrichment_individual.Rmd
- Pulls: counts crosstabs (syn64398049 and syn64398051); 
differential expression for MPNST, time, and drug pairs (syn69910880)
- Uploads: transcription factor activity for MPNST, time, and drug pairs (syn64420968)

## Figure 4A: Gene set correlation across all samples 
Table S3: 03_GSEA_individual.Rmd

- Pulls: differential expression for MPNST, time, and drug pairs (syn69910880)
- Uploads: enrichment of MSigDB hallmark gene sets for 
MPNST, time, and drug pairs (syn69910859)

## Figure 4B:
hand drawn network of targets and TFs

## Figure 4C: predicted target expression across three combinations

## Figure 4D DMEA prediction

## Drug mechanism enrichment analysis (DMEA)
Figure 4C, Table S4: 04_DMEA.Rmd
- Pulls: differential expression for time and drug pairs (syn64397743); 
public adherent CCLE RNA-seq and PRISM drug sensitivity
- Uploads: DMEA results based on sensitivity for time and drug pairs (syn65672258)

Table S4: 04_DMEA_individual.Rmd
- Pulls: differential expression for MPNST, time and drug pairs (syn69910880); 
public adherent CCLE RNA-seq and PRISM drug sensitivity
- Uploads: DMEA results based on sensitivity for MPNST, time and drug pairs (syn69910844)

### Evidence for MEK inhibitor combos
Figure 5: 05_MEKi_combos.R
- Pulls: PDX mirdametinib + vorinostat tumor size (syn68900596); 
differential expression (syn69910880); 
transcription factor activity (syn64693461); 
GSEA (syn69910859); DMEA (syn69910844) each for MPNST, time, and drug pairs

