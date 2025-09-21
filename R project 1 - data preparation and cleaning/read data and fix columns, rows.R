library(tidyverse)


setwd("C:/Users/willi/Desktop/C1 R code")
# Load gene data
gene_df <- read.table(file = 'data/count_matrix_genes.tsv', sep = '\t', header = TRUE)

head(gene_df)

# Fixing gene_df
ncol(gene_df)
colnames(gene_df)[2:length(gene_df)]

rownames(gene_df) <- gene_df$gene
gene_df <- gene_df[,2:ncol(gene_df)]


# =====================================================================================


setwd("C:/Users/willi/Desktop/C1 R code")
# Load TE data
TE_df <- read.table(file = 'data/count_matrix_transposable_elements.tsv', sep = '\t', header = TRUE)

head(TE_df)

# Fixing TE_df
ncol(TE_df)
colnames(TE_df)[2:length(TE_df)]

rownames(TE_df) <- TE_df$TE_ID
TE_df <- TE_df[,2:ncol(TE_df)]


setwd("C:/Users/willi/Desktop/C1 R code")
# Load participants
participant_data <-  read.table(file = 'data/participant_info_phenotype_data.tsv', sep = '\t', header = TRUE)
head(participant_data)

# Fixing participant_data
nrow(participant_data)
rownames(participant_data) <- participant_data$Individual
rownames(participant_data)
participant_data <- participant_data[,3:ncol(participant_data)]


colnames_TEdf <- substr(colnames(TE_df), 1, nchar(colnames(TE_df)) - 4)
colnames_gene_df <- substr(colnames(gene_df), 1, nchar(colnames(gene_df)))
rownames_participants <- rownames(participant_data)

common_samples <- intersect(colnames_TEdf, rownames_participants)
supercommon_samples <- intersect(colnames_gene_df, common_samples)

# Recode column names in TE_df
TE_df_columns_shortened <- TE_df
colnames(TE_df_columns_shortened) <- substr(colnames(TE_df_columns_shortened), 1, nchar(colnames(TE_df_columns_shortened)) - 4)

# Filter TE_df to keep only columns that are present in participant_data
TE_df_filtered <- TE_df_columns_shortened %>% select(any_of(supercommon_samples))
participant_data_filtered <- participant_data[rownames(participant_data) %in% supercommon_samples, ]
gene_df_filtered <- gene_df %>% select(any_of(supercommon_samples))

# Change all data types in matrix to integers:
TE_df_filtered <- TE_df_filtered %>%
  mutate(across(where(is.numeric), as.integer))

# Convert 'Status' (disease status) to a factor:
participant_data_filtered$Status <- factor(participant_data_filtered$Status, levels = c(1, 2))

rm(gene_df, participant_data, TE_df_columns_shortened, TE_df)



setwd("C:/Users/willi/Desktop/C1 R code/cleaned RData files")
saveRDS(file = 'gene_matrix.RDS', object = gene_df_filtered)
saveRDS(file = 'transposon_matrix.RDS', object = TE_df_filtered)
saveRDS(file = 'phenotype_data.RDS', object = participant_data_filtered)

print('ran to end of script')
