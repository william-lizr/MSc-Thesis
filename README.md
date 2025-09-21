## Contents of the file:

MSc-Thesis/
│
├── R project 1 - data preparation and/
│   ├── read data and fix columns, rows.R
│   └── splice_columns.R
│
├── R project 2 - descriptives and covariates/
│   ├── Ex 0 A PCA and identifying covariates.R
│   ├── Ex 0 B Descriptive stats dataset.R
│   ├── Ex 0 C NEW model comparisons.R
│   ├── Ex 0 C Testing for covariates format.R
│   └── Model comparison functions.R
│
├── R project 3 - base deseq model/
│   └── Differential expression.Rmd
│
├── R project 4 - gene enrichment analysis/
│   ├── gene enrichment script.Rmd
│   ├── glm count modelling rewrite.R
│   ├── numeric function.R
│   ├── test_1 -cpm-values.tsv
│   └── test_1 -lm-results.tsv
│
└── README.md


### Part I: Data preparation
Data read in from TSV and columns fixed for DESeq modelling

### Part II: Computing descriptives, inspecting covariates
A - Principal component analysis to assess certain covariates
B - Descriptive stats (mean and SD) for participant characteristics
C - Comparing whether models add signal based on LRT comparison within DESeq. Also comparing counts and log2fold change average, variance between models to see if covariates significantly change the model and should be included.
C - 'Testing for covariates'

### Part III:
Running differential expression using DESeq2 package using the chosen model.

### Part IV:
Gene enrichment based on proximity.
Step 1: Significantly differential expressed TEs (DE TEs) found using initial DESeq model. 
Step 2: Gene within 5000 bp of each TE found and paired with TEs in a one (TE) to many (genes) mapping.


