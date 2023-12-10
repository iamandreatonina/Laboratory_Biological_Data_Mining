# Set the working directory
setwd("D:/VarieTHOM/University/QCB/3_SEMESTRE/Data Mining/Laboratory (Blanzieri)/0_PROJECT/Datasets_finals")

# Load necessary libraries
library("enrichR")
library("ggplot2")
library("biomaRt")
library("clusterProfiler")
library("org.Hs.eg.db")

#### SKIP TO LOAD IF YOU HAVE LAREDY CONVERTED ####

# Read the dataset DEGs_subtype_T.csv and filter rows where 'class' is '+'
DEGs_subtype_T <- read.csv("DEGs_subtype_T.csv")
T_id <- DEGs_subtype_T[DEGs_subtype_T$class == "+", "X", drop = FALSE]

# Read the dataset DEGs_subtype_PT.csv and filter rows where 'class' is '+'
DEGs_subtype_PT <- read.csv("DEGs_subtype_PT.csv")
PT_id <- DEGs_subtype_PT[DEGs_subtype_PT$class == "+", "X", drop = FALSE]

# Read the dataset DEGs_subtype_B.csv and filter rows where 'class' is '+'
DEGs_subtype_PB <- read.csv("DEGs_subtype_B.csv")
PB_id <- DEGs_subtype_PB[DEGs_subtype_PB$class == "+", "X", drop = FALSE]

# Convert the ensembl to hugo

# Load the ensembl dataset for human genes
ensmebl <- useMart(biomart = 'ensembl', dataset = 'hsapiens_gene_ensembl')

# Convert ENSEMBL gene IDs to HUGO gene symbols for DEGs_subtype_T
convert_T <- getBM(attributes = c('ensembl_gene_id', 'entrezgene_id', 'external_gene_name'),
                   filters = c('ensembl_gene_id'), 
                   values = T_id$X, 
                   mart = ensmebl)

# Write the result to a CSV file
write.csv(convert_T, "convert_T.csv", row.names = TRUE)

# Convert ENSEMBL gene IDs to HUGO gene symbols for DEGs_subtype_PT
convert_PT <- getBM(attributes = c('ensembl_gene_id', 'entrezgene_id', 'external_gene_name'),
                    filters = c('ensembl_gene_id'), 
                    values = PT_id$X, 
                    mart = ensmebl)

# Write the result to a CSV file
write.csv(convert_PT, "convert_PT.csv", row.names = TRUE)

# Convert ENSEMBL gene IDs to HUGO gene symbols for DEGs_subtype_PB
convert_PB <- getBM(attributes = c('ensembl_gene_id', 'entrezgene_id', 'external_gene_name'),
                    filters = c('ensembl_gene_id'), 
                    values = PB_id$X, 
                    mart = ensmebl)

# Write the result to a CSV file
write.csv(convert_PB, "convert_PB.csv", row.names = TRUE)


#### Loead the file you alredy produced ####

convert_T <- read.csv("convert_T.csv")
convert_T$X<- NULL
convert_PT <- read.csv("convert_PT.csv")
convert_PT$X<- NULL
convert_PB <- read.csv("convert_PB.csv")
convert_PB$X<- NULL


# Connecting to Enrichr web service

# List available databases from Enrichr
dbs <- listEnrichrDbs()
dbs <- dbs[order(dbs$libraryName),]
Databases <- data.frame(dbs$libraryName)

# Enrichment analysis for DrugMatrix and IDG_Drug_Targets_2022 databases

# Define the databases for enrichment analysis
dbs_dd <- c("DrugMatrix", "IDG_Drug_Targets_2022")

# Perform enrichment analysis for DEGs_subtype_T
upClinical_T <- enrichr(genes = convert_T$external_gene_name, databases = dbs_dd)

# Perform enrichment analysis for DEGs_subtype_PT
upClinical_PT <- enrichr(genes = convert_PT$external_gene_name, databases = dbs_dd)

# Perform enrichment analysis for DEGs_subtype_PB
upClinical_PB <- enrichr(genes = convert_PB$external_gene_name, databases = dbs_dd)

# Extract results for DrugMatrix and IDG_Drug_Targets_2022 from DEGs_subtype_T ####
DM_up_T <- data.frame(upClinical_T[["DrugMatrix"]])
write.csv(DM_up_T, "DM_up_T.csv", row.names = TRUE)
# Sort the dataframe based on the "Combined.Score" column in descending order
sorted_DM_up_T <- head(DM_up_T[order(DM_up_T$Combined.Score, decreasing = TRUE), ],5)

IDG_up_T <- data.frame(upClinical_T[["IDG_Drug_Targets_2022"]])
write.csv(IDG_up_T, "IDG_up_T.csv", row.names = TRUE)
# Sort the dataframe based on the "Combined.Score" column in descending order
sorted_IDG_up_T <- head(IDG_up_T[order(IDG_up_T$Combined.Score, decreasing = TRUE), ],5)



#### Extract results for DrugMatrix and IDG_Drug_Targets_2022 from DEGs_subtype_PT ####
DM_up_PT <- data.frame(upClinical_PT[["DrugMatrix"]])
write.csv(DM_up_PT, "DM_up_PT.csv", row.names = TRUE)
# Sort the dataframe based on the "Combined.Score" column in descending order
sorted_DM_up_PT <- head(DM_up_PT[order(DM_up_PT$Combined.Score, decreasing = TRUE), ],5)

IDG_up_PT <- data.frame(upClinical_PT[["IDG_Drug_Targets_2022"]])
write.csv(IDG_up_PT, "IDG_up_PT.csv", row.names = TRUE)
### Sort the dataframe based on the "Combined.Score" column in descending order
sorted_IDG_up_PT <- head(IDG_up_PT[order(IDG_up_PT$Combined.Score, decreasing = TRUE), ],5)



# Extract results for DrugMatrix and IDG_Drug_Targets_2022 from DEGs_subtype_PB  ####
DM_up_PB <- data.frame(upClinical_PB[["DrugMatrix"]])
write.csv(DM_up_PB, "DM_up_PB.csv", row.names = TRUE)
# Sort the dataframe based on the "Combined.Score" column in descending order
sorted_DM_up_PB <- head(DM_up_PB[order(DM_up_PB$Combined.Score, decreasing = TRUE), ],5)

IDG_up_PB <- data.frame(upClinical_PB[["IDG_Drug_Targets_2022"]])
write.csv(IDG_up_PB, "IDG_up_PB.csv", row.names = TRUE)
# Sort the dataframe based on the "Combined.Score" column in descending order
sorted_IDG_up_PB <- head(IDG_up_PB[order(IDG_up_PB$Combined.Score, decreasing = TRUE), ],5)



plotEnrich(upClinical_T[[2]], showTerms = 20, numChar = 50, y = "Count", orderBy = "Combined.Score")

plotEnrich(upClinical_PT[[2]], showTerms = 20, numChar = 50, y = "Count", orderBy = "Combined.Score")

plotEnrich(upClinical_PB[[2]], showTerms = 20, numChar = 50, y = "Count", orderBy = "Combined.Score")




