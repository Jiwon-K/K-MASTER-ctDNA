# K-MASTER Pan-cancer ctDNA Analysis

This repository contains code used for genomic and statistical analyses in the following study:

> **Pan-cancer circulating tumor DNA atlas establishes a precision oncology framework for prognostic and therapeutic stratification**


The repository provides code for genomic alteration profiling, comparisons with tissue sequencing and an external ctDNA cohort, prognostic modeling, and treatment outcome analyses.


## Code description

- Figure1.R: K-MASTER ctDNA Genomic landscape, mutation count, MATH scores, and OncoKB annotations.
- Figure2.R: Comparison of genomic alteration frequencies between K-MASTER and MSK-ACCESS.
- Figure3.R: Comparison of genomic profiles between tissue and ctDNA.
- Figure4.R: Overall survival, Molecular Prognostic Index, elastic-net Cox regression, and SHAP analyses.
- Figure5.R: Time-to-next-treatment and pathway co-alteration analyses.


## Software

Analyses were performed using R version 4.3.3.

The analyses used R packages for genomic data processing, statistical modeling, survival analysis, and visualization. 


## Data availability

Patient-level clinical and genomic data from the K-MASTER cohort are deposited in the Clinical & Omics Data Archive (CODA; accession number CODA_D23008). 

Access to these data is controlled and requires approval through CODA. The controlled-access patient-level data are not included in this repository.


## Usage

The scripts require the corresponding input datasets. Data processing steps and required input files are described in the individual scripts.


## Contact:
For questions regarding the analysis code, please contact:

Jiwon Kim (greati4320@korea.ac.kr)
