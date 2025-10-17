# Project
MPNST-PDX-MT-drugViability

## Analysis scripts which read to and from Synapse
### Dose-response curve fitting and data compilation
Table S1: 00_createCurveStats.py
- Pulls: relative viability data (single agent: syn65473019 using metadata: syn65473034, combo: syn66330226 using metadata: syn66330284)
- Uploads: single agent relative viability table (syn65986622) and curve parameters (syn65941820); combo relative viability table (syn68156852)

### Single agent AUC heatmap
Figure 3A-B, 5A: 01_viabilityPlot.R
- Pulls: single agent viability curves (syn65941820)

### Synergy metric comparison and heatmaps
Figure S1B-C: 02_synergy.R
- Pulls: combo relative viability data (syn68156852); bliss synergy results (syn68639935); 
musyc synergy results (syn68736713); 48h median CI (syn70365485); 120h median CI (syn70365484); 
delta AUC (syn69978141)

### Dose-response curves
03_doseResponseCurves.R
- Pulls: single agent relative viability data (syn65986622); 
combo relative viability data (syn68156852); bliss synergy results (syn68900210); 
musyc synergy results (syn68736713)
- Uploads: compiled combo relative viability, bliss synergy, and musyc synergy results (syn68900322)

### Prediction comparison
04_predictionComparison.R
- Pulls: DMEA results (syn65672258) and related drug correlations (syn65672266); compiled combo relative viability, bliss synergy, and musyc synergy results (syn68900322)
- Uploads: shared DMEA and synergy results (syn68900547); shared drugs across predictions and synergy scores (syn68900548)

### PDX studies
Figure 5B: 05_PDX.R
- Pulls: mean PDX tumor sizes with mirdametinib, vorinostat, mirdametinib+vorinostat, or vehicle (MN-2: syn69953578, WU-225: syn68900596)
- Saves: plots of mean tumor size over time for each PDX

