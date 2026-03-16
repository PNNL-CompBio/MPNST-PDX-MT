# Leveraging a Precision Medicine Platform to Predict Novel Therapies for Malignant Peripheral Nerve Sheath Tumors

This project contains analysis for the DoD funded project across University of Minnesota, Johns Hopkins University, 
Washington University of St. Louis, and Pacific Northwest National Laboratory to develop a platform for idnetification of drug combinations 
in NF1 MPNSTs.

## The EXACT Pipeline
Currently we have built a pipeline that leverages MPNST PDX in an 3D culture system to characterize
drug response to enable selection of drug combinations. For more information please see the [EXACT project directory](./EXACT-pipeline). 

## Data Types
Here we have scripts required to process the three primary data modalities we are collecting including:
1. **Drug viability data:** Here we have dose/response values that we use to assess Area Under the Curve (AUC) and
other drug metrics for both single and combination screens in our 3D Meds Platform. This analysis has now been folded into our 
EXACT platform and can be found in [./EXACT-pipeline/drugViability](./EXACT-pipeline/drugViability) directory.
2. **RNA Sequencing: ** We have collected dozens of responses to various drugs in our 3D Meds platform. Initial processing scripts can be found in
the [RNASeq](./RNASeq) directory and the final processing can be found in the [./EXACT-pipeline/RNASeq](./EXACT-pipeline/RNASeq) 
3. **Proteomics: ** We are in the process of optimizing our proteomic analysis for drug response measurement. Preliminary analysis can be found 
in the [proteomics](./proteomics/) directory. 


## Project Data
All data for this project is stored on Synapse at
http://synapse.org/mpnst. To run the analysis described here you will need to request data access.

