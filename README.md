# Sepsis-MR

This repository contains the data files, R code and figures for the paper "Cardiometabolic traits, sepsis and severe covid-19: a Mendelian randomization investigation" (M. Ponsford et al. 2020). A link to the paper will appear here soon.

We conduct Mendelian randomization to assess the existence of a causal relationship between five cardiometabolic traits and the risk of sepsis, as well as with severe covid-19 infection risk. Specifically, we consider the following traits:
- Body mass index
- LDL cholesterol
- Systolic blood pressure
- Smoking
- Type 2 diabetes mellitus

Genetic associations with each trait were obtained from relevant GWAS consortia - see our manuscript for references. They are included here in the excel files named after each trait.

Genetic associations with sepsis risk were obtained from the UK Biobank and the HUNT longitudinal study. Summary-level data are included in the text files "Sepsis Data from UKBB" and "Sepsis Data from HUNT" respectively. Genetic associations with severe covid-19 infection risk were obtained from a study recently oublished at NEJM - see our manuscript for details. Summary-level data are included in the text file "Covid Data" . 

The R file "Sepsis MR" contains code to read the data into R, conduct Mendelian randomization analysis, report results and create the three figures included in the paper. The MR analysis is conducted using the R package "MendelianRandomization" and is run separately for UK Biobank and HUNT data. Likewise, the R file "Covid MR" contains code to import, analyse and plot the covid-19 data.

The figures included in the manuscript have also been uploaded in this repository for transparency.

