# Project
MPNST-PDX-MT-drugViability

## Analysis scripts which read to and from Synapse
### Dose-response curve fitting and data compilation
00_createCurveStats.py
- Pulls: relative viability data (single agent: syn65473019 using metadata: syn65473034, combo: syn66330226 using metadata: syn66330284)
- Uploads: single agent relative viability table (syn65986622) and curve parameters (syn65941820); combo relative viability table (syn68156852)

### Single agent AUC heatmap
01_viabilityPlot.R
- Pulls: single agent viability curves (syn65941820)

### Synergy metric comparison and heatmaps
02_synergy.Rmd
- Pulls: combo relative viability data (syn68156852); bliss synergy results (syn68639935); musyc synergy results (syn68736713)

### Dose-response curves
03_doseResponseCurves.Rmd
- Pulls: single agent relative viability data (); combo relative viability data (syn68156852); bliss synergy results (syn68900210); musyc synergy results (syn68736713)
- Uploads: compiled combo relative viability, bliss synergy, and musyc synergy results (syn68900322)

### Prediction comparison
04_predictionComparison.R
- Pulls: DMEA results (syn65672258) and related drug correlations (syn65672266); compiled combo relative viability, bliss synergy, and musyc synergy results (syn68900322)
- Uploads: shared DMEA and synergy results (syn68900547); shared drugs across predictions and synergy scores (syn68900548)

