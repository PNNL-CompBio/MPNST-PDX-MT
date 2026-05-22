## Transcript-drug response analysis

This directory studies various ways to compare the gene expression values
and the drug response as well as the combination drug response. 

There is a single helper file, `correlation_functions.R`,
that contains functions to carry out the markdowns listed below. We employed 
two distinct approaches to drug response: one that is correlation based
and one that is regression based. For each of those, we used both the
genes as features and ssGSEA scores (using cancer hallmarks)

## AUC and AUC shift
The first thing we do to compare treated gene expression to the drug response
is to determine a metric of drug response in combination. For this we define 
the AUC shift to be the `AUC_d - AUC_initial`. Here we now compare the effect
a specific drug `d` has on a combination across initial drugs.

[`00_aucExploration.Rmd`](./00_aucExploration.Rmd) plots the AUC and AUC shift
of individual drugs in our samples.

## Basic correlation analysis
Here we used the 'AUC shift' for each drug and asked how much the baseline 
expression differences correlated with it to ask if individual features could 
give us a sense of downstream drug effect.

### Gene based
 
This is explored in 
[`01_geneCorrelationEnrichment.Rmd`](./01_geneCorrelationEnrichment.Rmd). 
The issue with this ability is that the pathway
enrichment comes after the correlation, and there are many overlapping pathways.

### Pathway based
To identify which pathways were uniquely correlated with drug response we
used ssGSEA scores on the samples themselves to ask which pathways were
correlated with AUC.


The file [`02_gseaCorrelation.Rmd`](./02_gseaCorrelation.Rmd) explores this. This
file is still under development as we work to identify more trends.

## Regression and prediction

We used the single agent drug response values to build a very small elastic
net model of drug response for the drugs of interest, then applied it to the 
treated samples.

_TODO_: We need to augment this model using additional training data!

### Gene based regression
Here we used an elastic net model to regress against all genes for a single
AUC as well as the AUC shift. The file 
[`03_geneRegressionTests.Rmd`](./03_geneRegressionTets.Rmd) shows our early results. 

### Pathway based regression

Here we used an elastic net model to regress against all pathways for a single
AUC as well as the AUC shift. The file 
[`04_gseaRegressionTests.Rmd`](./03_geneRegressionTets.Rmd) shows our early results. 
