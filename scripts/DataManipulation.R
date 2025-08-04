#loading library
library(tidyverse)
library(GEOquery)
library(openxlsx)
library(ggpubr)
library(DESeq2)

#get Metadata
data <- getGEO(GEO = "GSE130078")
metadata <- pData(phenoData(data[[1]]))
metadata_subset <- metadata |> 
  select(c(1,2,39))
metadata_modified <- metadata_subset |> 
  rename("sample_type:ch1" = "tissue") |> 
  mutate(title = gsub("_[^_]+$","",title))
metadata_modified$tissue <- factor(metadata_modified$tissue, levels = c("normal", "cancer"))

#Creating a unified count matrix
#1. getting a list of all the text file's addresses
count_files <- list.files(path = "data", pattern = "\\.assigned_count\\.txt$", full.names = TRUE)

#2. converting each text files into data frames
count_list <- lapply(count_files, function(file){
  read.delim(file, header = TRUE, stringsAsFactors = FALSE)
})

#3. creating a vector for each data frame
count_named_list <- lapply(count_list,function(df){
  count_number <- as.numeric(df$Assigned_count)
  names(count_number) <- df$Gene_id
  count_number
})

#4. Finally creating the unified count matrix
count_matrix <- do.call(cbind, count_named_list)
colnames(count_matrix) <- c(metadata$geo_accession)
count_df <- as.data.frame(count_matrix)
colnames(count_df) <- c(metadata$geo_accession)

#Normalizing the count
dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = metadata_modified, design = ~ tissue)
dds <- estimateSizeFactors(dds)
head(dds)
normalized_counts<-counts(dds, normalized=TRUE)
normalized_counts_df <- as.data.frame(normalized_counts)

#creating a long sample
normalized_counts_df$gene_id <- rownames(normalized_counts_df)
counts_data_long <- normalized_counts_df|> pivot_longer(-gene_id, names_to = "Samples", values_to = "fpkm")

#joining meta data and long sample
combined_matrix <- counts_data_long |> left_join(metadata_modified, by = c("Samples" = "geo_accession"))
write.csv(combined_matrix, "data/GSE130078_counts.csv", row.names = FALSE)

