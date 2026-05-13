# Ex vivo assessment of combination therapies (EXACT): An integrated platform to discover combination cancer treatments
The goal of this project is to leverage a diverse set
of 3D-MEDS models grown from patient-derived xenografts (PDX)
to identify potential drug combinations that could be used to treat
MPNST patients.Specifically we leverage dose response data and gene expression
measurements from these models to prioritize drug combinations.

## Project Data
All data for this project is stored on Synapse at
http://synapse.org/mpnst. It is currently under embargo as we continue
to process and analyze the data but will release the data upon
publication. 

## Analysis scripts
Below are the basic steps to analyze the data and reproduce the figures in the 
Shah et al manuscript.

### Install dependencies
Dependencies can be installed by running installDependencies.R

### Single and combo drug measurements
The first half of the EXACT analysis is derived from the drug viability scores. 
These are calculated and plotted in the [drugViability](drugViability/) directory.
This also includes code for Figures 2A and 2B.

### Single drug transcriptional response
These data are currently being analyzed in our [RNAseq](RNAseq/)
directory. 

### Supplementary table compilation
Supplementary table 2: suppTable2Compilation.R
