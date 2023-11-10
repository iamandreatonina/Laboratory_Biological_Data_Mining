library(dplyr)

# Specify the directory
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
write.csv(final_data2, file = "merged_GSM181157.csv", row.names = FALSE)
