## Contents of the file:

MSc-Thesis/  <br> 
│  <br> 
├── R project 1 - data preparation and cleaning/  <br> 
│   ├── read data and fix columns, rows.R  <br> 
│   └── splice_columns.R  <br> 
│<br> 
├── R project 2 - descriptives and covariates/<br> 
│   ├── Ex 0 A PCA and identifying covariates.R <br> 
│   ├── Ex 0 B Descriptive stats dataset.R <br> 
│   ├── Ex 0 C NEW model comparisons.R <br> 
│   ├── Ex 0 C Testing for covariates format.R <br> 
│   └── Model comparison functions.R <br> 
│ <br> 
├── R project 3 - base deseq model/ <br> 
│   └── Differential expression.Rmd <br> 
│ <br> 
├── R project 4 - gene enrichment analysis/ <br> 
│   ├── gene enrichment script.Rmd <br>  
│   ├── glm count modelling rewrite.R <br> 
│   ├── numeric function.R <br> 
│   ├── test_1 -cpm-values.tsv <br> 
│   └── test_1 -lm-results.tsv <br> 
│ <br> 
└── README.md <br> 


#### Part I: Data preparation <br> 
Data read in from TSV and columns fixed for DESeq modelling <br> 
 
### Part II: Computing descriptives, inspecting covariates <br> 
##### A - Principal component analysis to assess certain covariates <br> 
##### B - Descriptive stats (mean and SD) for participant characteristics <br> 
##### C - Comparing whether models add signal based on LRT comparison within DESeq. Also comparing counts and log2fold change average, variance between models to see if covariates significantly change the model and should be included. <br> 
##### C - 'Testing for covariates' <br> 

#### Part III:
Running differential expression using DESeq2 package using the chosen model.

#### Part IV:
Gene enrichment based on proximity. <br> 
##### Step 1: Significantly differential expressed TEs (DE TEs) found using initial DESeq model.  <br> 
  Data object: TE characteristics (family, subfamily, locus) | genomic coordinate start | genomic coordinate end | chromosome | significance | adjusted p-value (filtered to TEs < 0.05 adjusted alpha) <br> 
##### Step 2: Genes within 5000 bp of each TE found and paired with TEs in a one (TE) to many (genes) mapping. <br> 
  Data object: Significant TE (family | subfamily | locus | coordinates | chromosome) | Gene within 5000 bp (ENCODE ID | genomic coordinate start | genomic coordinate end )
  Number of rows: for each TE find N genes within 5000bp <br> 
##### Step 2a: Filter expression matrix
Filter expression matrix to only include expression data of genes which are relevant to the analysis.
##### Step 3: Correlation analysis of each-TE gene pair
Assess significance of correlation between TE expression (across healthy and ALS patients) and gene expression (across healthy and ALS patients).
Summary of counts of genes which are significantly correlated after multiple comparisons yields information about which pathways may be altered in ALS.

Binomial distribution (like DESeq2) used for correlation analysis.
Multiple comparisons correction done using Benjamini–Hochberg procedure

  


