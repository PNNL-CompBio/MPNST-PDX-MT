# Leveraging a Precision Medicine Platform to Predict Novel Therapies for Malignant Peripheral Nerve Sheath Tumors

This project contains analysis for the DoD funded project across University of Minnesota, Johns Hopkins University, 
Washington University of St. Louis, and Pacific Northwest National Laboratory to develop a platform for identification of drug combinations 
in NF1 MPNSTs.

This repository includes both basic data processing steps and also the code included as part of our EXACT project.

## Data Processing
Here are scripts that are used to process and explore the data:
1. **RNA Sequencing: ** We have collected dozens of responses to various drugs in our 3D Meds platform. Initial processing scripts can be found in
the [RNASeq](./RNASeq) directory and the final processing can be found in the [./EXACT-pipeline/RNASeq](./EXACT-pipeline/RNASeq) 
3. **Proteomics: ** We are in the process of optimizing our proteomic analysis for drug response measurement. Preliminary analysis can be found 
in the [proteomics](./proteomics/) directory. 


## The EXACT Pipeline
The pathway and dose response curve analysis that is part of the Shah et al
manuscript to integrate dose response curve data with gene expression data
to identify putative drug combinations. This is described in the 
[EXACT-pipeline](./EXACT-pipeline) directory.

## Project Data

All data for this project is stored on Synapse at
http://synapse.org/mpnst. To run the analysis described here you will need to request data access.

