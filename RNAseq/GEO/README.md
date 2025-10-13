## Analysis of MEK inhibitor + HDAC inhibitor transcriptomic 
## datasets in PRC2-null MPNSTs found on Gene Expression Omnibus (GEO)
# Differential expression analyses
00_MEKi-treated-MPNST-diffexp-24h.R
- Pulls: counts tables from GEO (GSE183307 and GSE262030)
- Saves: differential expression (GSE183307_diffexp.csv and GSE262030_diffexp.csv)

# Drug mechanism enrichment analysis (DMEA): sensitivity predictions
01_DMEA_sensitivity.R
- Pulls: differential expression (GSE183307_diffexp.csv and GSE262030_diffexp.csv)
- Saves: drug sensitivity predictions (DMEA_results.csv and Compiled_DMEA_results.csv) and plots

# Drug mechanism enrichment analysis (DMEA): similarity predictions
02_DMEA_similarity.R
- Pulls: differential expression (GSE183307_diffexp.csv and GSE262030_diffexp.csv)
- Saves: drug similarity predictions (DMEA_Cmap_L1000.csv) and plots; also saves 
drug and mechanism of action (MOA) rankings produced by Cmap L1000 Query tool 
(Cmap_L1000_drugRanks.csv and Cmap_L1000_moaRanks.csv) but I use DMEA_Cmap_L1000.csv 
for consistent calculations with sensitivity predictions

# Result compilation: drug sensitivity and similarity predictions
03_resultCompilation.R
- Pulls: drug sensitivity predictions (DMEA_results.csv) and 
drug similarity predictions (DMEA_Cmap_L1000.csv)
- Saves: plots (venn diagrams and dot plots)
