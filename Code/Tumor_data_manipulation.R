library(dplyr)
library(EDASeq)
library(biomaRt)
library(edgeR)
library(readxl)

##### Make a unique merged dataframe for GSE181157

# Specify the directory path 
directory_path <- "D:/VarieTHOM/University/QCB/3_SEMESTRE/Data Mining/Laboratory (Blanzieri)/0_PROJECT/GSE181157_RAW"

# List all files that end with ".txt"
file_names <- list.files(directory_path, pattern = "\\.txt$", full.names = TRUE)

# Extract file names without extension and remove ".counts" part
file_names_cleaned <- sub("\\.counts$", "", tools::file_path_sans_ext(basename(file_names)))

# Initialize
final_data2 <- data.frame()

# Loop through each file
for (i in seq_along(file_names)) {
  # Read the file into a temporary data frame
  temp_data <- read.table(file_names[i], col.names = c("ENSEMBL_ID", file_names_cleaned[i]))
  
  # Merge the temporary data frame with the final data frame
  if (nrow(final_data2) == 0) {
    final_data2 <- temp_data
  } else {
    final_data2 <- merge(final_data2, temp_data, by = "ENSEMBL_ID", all = TRUE)
  }
}

# View the final data frame
View(final_data2)

# Write the final data frame to a text file
write.csv(final_data2, file = "merged_GSE181157.csv", row.names = FALSE)



##### Data manipulation GSE227832
GSE227832 <- read.csv('GSE227832_RNAseq_read_counts.txt',header = T,sep = '\t')

# Now we select just the useful data by eliminating the others 
GSE227832 <- GSE227832 %>% dplyr::select(-c(332:341)) %>% 
              dplyr::select(-c(ALL_validation_12,ALL_Validation_20,Constitutional_458rep3,Constitutional_559rep2,ALL_317_2,ALL_317_3,ALL_468_2,ALL_555_2,ALL_680_2))

# substitution of first column with ensembl_gene_id
colnames(GSE227832)[1] <- 'ensembl_gene_id'

#####  Data manipulation GSE133499

GSE133499 <- readxl::read_xlsx('GSE133499_count_matrix_samples_sciarrillo.xlsx')

# remove useless data 

GSE133499 <- GSE133499 %>% dplyr::select(-c(40:43))

# create the ensembl object that points to the Hsapiens database

ensembl_GSE133499 <- biomaRt::useEnsembl(biomart = "ensembl",dataset = "hsapiens_gene_ensembl")



