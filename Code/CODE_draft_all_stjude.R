setwd("D:/VarieTHOM/University/QCB/3_SEMESTRE/Data Mining/Laboratory (Blanzieri)/0_PROJECT/DATAsets")

# Read the data from text files
data_mutations_2016 <- read.table("all_stjude_2016/data_mutations.txt", header = TRUE, sep = "\t", fill = TRUE)
meta_mutation_2016 <- read.table("all_stjude_2016/meta_mutations.txt", header = TRUE, sep = "\t")
data_clinical_sample_2016 <- read.table("all_stjude_2016/data_clinical_sample.txt", header = TRUE, sep = "\t")
data_clinical_patient_2016 <- read.table("all_stjude_2016/data_clinical_patient.txt", header = TRUE, sep = "\t")

# Load the dplyr package
library(dplyr)

# Use dplyr to drop columns with all NA values in the data_mutations data frame
data_filtered_2016 <- data_mutations_2016 %>%
  select_if(~!all(is.na(.)))

# Display the first few rows of the filtered data frame for inspection
tibble(data_filtered_2016)

# Save the data_filtered data frame in CSV format
write.csv(data_filtered_2016, "all_stjude_2016/data_filtered_2016.csv", row.names = FALSE)

library(readr)
filtered_2016 <- read_csv("all_stjude_2016/data_filtered_2016.csv")
tibble(filtered_2016)
###############################################################################################################

# Read the data from text files
data_mutations_2015 <- read.table("all_stjude_2015/data_mutations.txt", header = TRUE, sep = "\t", fill = TRUE)
meta_mutation_2015 <- read.table("all_stjude_2015/meta_mutations.txt", header = TRUE, sep = "\t")
data_clinical_sample_2015 <- read.table("all_stjude_2015/data_clinical_sample.txt", header = TRUE, sep = "\t")
data_clinical_patient_2015 <- read.table("all_stjude_2015/data_clinical_patient.txt", header = TRUE, sep = "\t")


# Use dplyr to drop columns with all NA values in the data_mutations data frame
data_filtered_2015 <- data_mutations_2015 %>%
  select_if(~!all(is.na(.)))

# Display the first few rows of the filtered data frame for inspection
tibble(data_filtered_2015)

# Save the data_filtered data frame in CSV format
write.csv(data_filtered_2015, "all_stjude_2015/data_filtered_2015.csv", row.names = FALSE)

library(readr)
filtered_2015 <- read_csv("all_stjude_2015/data_filtered_2015.csv")
tibble(filtered_2015)
