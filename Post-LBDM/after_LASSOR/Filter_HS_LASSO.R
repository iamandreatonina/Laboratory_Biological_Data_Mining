# Read data from "HS_after_LASSO.txt" into a DataFrame
data <- read.table("HS_after_LASSO.txt", header = TRUE, sep = "\t")

# Read data from "Origin_HS.csv" into another DataFrame
HS <- read.csv("Origin_Pre_B_HS.csv")

# Check if values in data$HS are present in Origin_HS$Ensembl.ID
values_in_Origin_HS <- data$HS %in% HS$Ensembl.ID

# Print the result (a logical vector indicating presence or absence)
print(values_in_Origin_HS)

# Filter Origin_HS based on values in data$HS
filtered_Origin_HS <- HS[HS$Ensembl.ID %in% data$HS, ]

# Write the filtered DataFrame to a new CSV file
write.csv(filtered_Origin_HS, file = 'LASSO_Origin_Pre_B_HS.csv', row.names = TRUE)

rm(data,HS,values_in_Origin_HS,filtered_Origin_HS)

