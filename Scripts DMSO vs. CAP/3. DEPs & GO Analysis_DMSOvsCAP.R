install.packages("dplyr") 
library(dplyr)

# Upload the .csv file
report.pg_matrix <- read.csv ("DE_results_DMSOvsCAP.csv",
                              sep = ";", 
                              header = TRUE, 
                              dec = ".", 
                              check.names = FALSE) 
# Filtering differentially expressed proteins
filtered_report.pg <- report.pg_matrix %>%
  filter(abs(log2_FC) > 1 & p_val < 0.05)
# Save the filtered file as a new CSV file
write.csv(filtered_report.pg, "DE_results_filtered.csv", row.names = FALSE)

# Getting protein IDs for further enrichment of pathways
DE_results_filtrado <- read.csv ("DE_results_filtered.csv",
                              sep = ";", 
                              header = TRUE, 
                              dec = ".", 
                              check.names = FALSE) 
protein_ids <- DE_results_filtrado$Protein_ID

# Install package ‘clusterProfiler’.
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler", force = TRUE)
library(clusterProfiler)

# Install package org.Hs.eg.db
BiocManager::install("org.Hs.eg.db", force = TRUE)
library(org.Hs.eg.db)  # Database for humans

# Convert genes to Entrez ID
entrez_ids <- bitr(protein_ids, fromType = "UNIPROT", 
                   toType = "ENTREZID", 
                   OrgDb = org.Hs.eg.db)

# Remove duplicates, keeping only the first mapping
entrez_ids_unique <- entrez_ids[!duplicated(entrez_ids$UNIPROT), ]

# Enrichment of GO pathways

# Install necessary packages
install.packages("BiocManager")   # If you do not have Bioconductor installed
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot", "ggplot2", force = TRUE))
library(clusterProfiler)  # For enrichment analysis
library(org.Hs.eg.db)     # Human gene annotation database
library(enrichplot)       # For display
library(ggplot2)          # Additional graphics

# Gene list in ENTREZ format
genes_entrez <- entrez_ids$ENTREZID  

# Conduct GO (Biological Process) enrichment analysis
ego_BP <- enrichGO(gene = genes_entrez, 
                   OrgDb = org.Hs.eg.db, 
                   keyType = "ENTREZID",
                   ont = "BP",  # "BP" para Biological Process, "MF" para Molecular Function, "CC" para Cellular Component
                   pAdjustMethod = "BH", 
                   pvalueCutoff = 0.05, 
                   qvalueCutoff = 0.05, 
                   readable = TRUE)  # Converts ENTREZ to readable gene names

# Dot plot 
dotplot(ego_BP, showCategory = 15) + 
ggtitle("Top 15 GO Biological Process terms")

# Bar plot
barplot(ego_BP, showCategory = 15, title = "GO Enrichment Analysis")

# Network diagram (to see relationships between enriched terms)
cnetplot(ego_BP, categorySize = "pvalue", foldChange = NULL)
