setwd("D:/VarieTHOM/University/QCB/3_SEMESTRE/Data Mining/Laboratory (Blanzieri)/0_PROJECT/Datasets_finals/ML_HS")
library("dplyr")

#### Read and process each CSV file ####

Predict_Knn <- read.csv("Predict_Knn.csv", header  = TRUE)
row.names(Predict_Knn) <- Predict_Knn$X
Predict_Knn$X <- NULL

Predict_nb <- read.csv("Predict_nb.csv", header  = TRUE)
row.names(Predict_nb) <- Predict_nb$X
Predict_nb$X <- NULL

Predict_rf <- read.csv("Predict_rf.csv", header  = TRUE)
row.names(Predict_rf) <- Predict_rf$X
Predict_rf$X <- NULL

Predict_xgb <- read.csv("Predict_xgb.csv", header  = TRUE)
row.names(Predict_xgb) <- Predict_xgb$X
Predict_xgb$X <- NULL


# Sample list of names
name_list1 <- my_levels$type
name_list2 <- my_levels$C_T
name_list3 <- my_levels$Cell_type
name_list4 <- my_levels$Cell_type[1:4]


#### Predicted XGB ####

Predict_xgb <- Predict_xgb %>%
  mutate(type = case_when(
    type == 1 ~ name_list1[1],
    type == 2 ~ name_list1[2]
  ))


Predict_xgb <- Predict_xgb %>%
  mutate(C_T = case_when(
    C_T == 1 ~ name_list2[1],
    C_T == 2 ~ name_list2[2]
  ))


Predict_xgb <- Predict_xgb %>%
  mutate(Cell_type = case_when(
    Cell_type == 1 ~ name_list3[1],
    Cell_type == 2 ~ name_list3[2],
    Cell_type == 3 ~ name_list3[3],
    Cell_type == 4 ~ name_list3[4],
    Cell_type == 5 ~ name_list3[5],  
  ))


Predict_xgb <- Predict_xgb %>%
  mutate(predicted_xgb = case_when(
    predicted_xgb == 1 ~ name_list4[1],
    predicted_xgb == 2 ~ name_list4[2],
    predicted_xgb == 3 ~ name_list4[3],
    predicted_xgb == 4 ~ name_list4[4] 
  ))




#### Predicted KNN ####

Predict_Knn <- Predict_Knn %>%
  mutate(type = case_when(
    type == 1 ~ name_list1[1],
    type == 2 ~ name_list1[2]
  ))


Predict_Knn <- Predict_Knn %>%
  mutate(C_T = case_when(
    C_T == 1 ~ name_list2[1],
    C_T == 2 ~ name_list2[2]
  ))


Predict_Knn <- Predict_Knn %>%
  mutate(Cell_type = case_when(
    Cell_type == 1 ~ name_list3[1],
    Cell_type == 2 ~ name_list3[2],
    Cell_type == 3 ~ name_list3[3],
    Cell_type == 4 ~ name_list3[4],
    Cell_type == 5 ~ name_list3[5],  
  ))


Predict_Knn <- Predict_Knn %>%
  mutate(predicted_Knn = case_when(
    predicted_Knn == 1 ~ name_list4[1],
    predicted_Knn == 2 ~ name_list4[2],
    predicted_Knn == 3 ~ name_list4[3],
    predicted_Knn == 4 ~ name_list4[4] 
  ))





#### unify the results ####
# Create a new data frame 'All_results' by selecting only the 'Cell_type' column from 'my_data'
All_results <- subset(my_data, select = Cell_type)  

# Add predicted values from the Random Forest model to 'predicted_RF' column
All_results$predicted_RF <- Predict_rf$predicted_rf

# Add predicted values from the K-Nearest Neighbors model to 'predicted_KNN' column
All_results$predicted_KNN <- Predict_Knn$predicted_Knn

# Add predicted values from the XGBoost model to 'predicted_XGB' column
All_results$predicted_XGB <- Predict_xgb$predicted_xgb

# Add predicted values from the Naive Bayes model to 'predicted_NB' column
All_results$predicted_NB <- Predict_nb$predicted_nb



#### Encode and compare with basic functions ####

columns_to_encode <- c('Cell_type','predicted_RF', 'predicted_KNN', 'predicted_XGB', 'predicted_NB')

# Convert specified columns to factor type
All_results <- All_results %>% mutate_at(columns_to_encode, as.factor)
setwd("D:/VarieTHOM/University/QCB/3_SEMESTRE/Data Mining/Laboratory (Blanzieri)/0_PROJECT/Datasets_finals")
write.csv(All_results, "All_results.csv", row.names = TRUE)

summary(All_results[, c('Cell_type','predicted_RF', 'predicted_KNN', 'predicted_XGB', 'predicted_NB')])

pairs(All_results[, c('Cell_type','predicted_RF', 'predicted_KNN', 'predicted_XGB', 'predicted_NB')])


#### Compare with a logical column ####

# Selecting the relevant columns
prediction_columns <- c('predicted_RF', 'predicted_KNN', 'predicted_NB')

# Check if values are the same across columns for each row
all_results_same_values <- apply(All_results[, prediction_columns], 1, function(row) all(row == row[1]))

# Create a logical column indicating if all predictions are the same for each row
All_results$all_predictions_same <- all_results_same_values


# create a df with only the models of interest and write it as csv
Compared_RF_KNN_NB <- subset(All_results, select = c('Cell_type','predicted_RF', 'predicted_KNN','predicted_NB','all_predictions_same'))  
setwd("D:/VarieTHOM/University/QCB/3_SEMESTRE/Data Mining/Laboratory (Blanzieri)/0_PROJECT/Datasets_finals/ML_HS")
write.csv(Compared_RF_KNN_NB, "Compared_RF_KNN_NB.csv", row.names = TRUE)


# Calculate the percentage of rows where all predictions are the same
percentage_same_values <- mean(all_results_same_values) * 100


# Create a color palette
colors <- c("red", "green")

# Create a Plotly bar chart
plot_ly(
  type = "bar",
  x = c("Same Predictions", "Different Predictions"),
  y = c(percentage_same_values, 100 - percentage_same_values),
  marker = list(color = colors)
) %>%
  layout(
    title = "Percentage of Rows with Same Predictions",
    yaxis = list(title = "Percentage", range = c(0, 100)),
    xaxis = list(title = "Prediction Type"),
    barmode = "stack"
  )

#pairs(All_results[, c('Cell_type','predicted_RF', 'predicted_KNN', 'predicted_XGB', 'predicted_NB')])
summary(All_results[, c('Cell_type','predicted_RF', 'predicted_KNN', 'predicted_XGB', 'predicted_NB')])
summary(Compared_RF_KNN_NB[,c('predicted_RF', 'predicted_KNN', 'predicted_NB')])






