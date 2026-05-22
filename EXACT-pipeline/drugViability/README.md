# EXACT Drug viability

### Dose-response curve fitting and data compilation
The first script [00_createCurveStats.py] collects the dose response measurements from Synapse and 
calculates the dose response data, which is uploaded as Table S1. 
- Pulls: relative viability data (single agent: syn65473019 using metadata: syn65473034, 
combo: syn66330226 using metadata: syn66330284)
- Uploads: single agent relative viability table (syn65986622) and curve 
parameters (syn65941820); combo relative viability table (syn68156852)

### Single agent AUC heatmap
[ 01_viabilityPlot.R] calculates the AUC and viability of the single agent drugs 
and generates Figure 3A-B.
- Pulls: single agent viability curves (syn65941820)

### Synergy metric comparison and heatmaps
[02_synergy.R] compares synergy results across different metrics and generates Figure S1: 
- Pulls: combo relative viability data (syn68156852); bliss synergy results (syn68639935); 
musyc synergy results (syn68736713); 48h median CI (syn70365485); 120h median CI (syn70365484); 
delta AUC (syn69978141)

### Dose-response curves
We also included code to plot dose response curves in [03_doseResponseCurves.R]

Creates the double agent dose response curves needed for Figure 7A.
- Pulls: single agent relative viability data (syn65986622); 
combo relative viability data (syn68156852); bliss synergy results (syn68900210); 
musyc synergy results (syn68736713)
- Uploads: compiled combo relative viability, bliss synergy, and musyc synergy results (syn68900322)

### PDX studies
Here we plot the curve for mirdametinib + vorinostat drug combination in two PDX
in [04_PDX.R].
This creates Figure 7B.
- Pulls: mean PDX tumor sizes with mirdametinib, vorinostat, mirdametinib+vorinostat, or vehicle (MN-2: syn69953578, WU-225: syn68900596)
- Saves: plots of mean tumor size over time for each PDX

