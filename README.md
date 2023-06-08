# Immune microenvironment profiling of normal appearing colorectal mucosa biopsied over repeat patient visits reproducibly separates lynch syndrome patients based on their history of colon cancer

This repository contains the code and data for our FrontiersIn manuscript as titled above

The repositary includes two .R scripts: **rxCOV.R** and **biomarkerSelection.R**. The rxCOV.R script identifies which cytokines have high signal fidelity for inclusion in further analysis. The biomarkerSelection.R takes the high fidelity cytokines and selects those ones that are best able to separate the LS patient groups being compared.

It also includes .RData files: rxCOVData.RData and patientData.RData. These data files contain the data used for the analysis and results presented in the manuscript.

Data file rxCOVData.RData is the input to the rxCOV.R script; and data file patientData.RData is the input to the biomarkerSelection.R script.


