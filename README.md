# Sepsis-MR

This repository contains the data files, R code and figures for the paper "Cardiometabolic traits and risk of sepsis: a Mendelian randomization study" (M. Ponsford et al. 2020). A link to the paper will appear here in due time.

We conduct Mendelian randomization to assess the existence of a causal relationship between five cardiometabolic traits and the risk of sepsis. Specifically, we consider the following traits:
- Body mass index
- LDL cholesterol
- Systolic blood pressure
- Smoking
- Type 2 diabetes mellitus

Genetic associations with each trait were obtained from relevant GWAS consortia - see the manuscript for references. They are included here in the excel files named after each trait.

Genetic associations with sepsis risk were obtained from the UK Biobank and the HUNT longitudinal study. Summary-level data are included in the text files "Sepsis Data from UKBB" and "Sepsis Data from HUNT" respectively.

The R file "Sepsis MR" contains code to read the data into R, conduct Mendelian randomization analysis, report results and create the three figures included in the paper. The MR analysis is conducted using the R package "MendelianRandomization" and is run separately for UK Biobank and HUNT data.

The three figures included in the manuscript have also been uploaded in this repository for transparency.

