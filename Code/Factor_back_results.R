setwd("D:/VarieTHOM/University/QCB/3_SEMESTRE/Data Mining/Laboratory (Blanzieri)/0_PROJECT/Datasets_finals")
library("dplyr")

#### Read and process each CSV file ####

Predict_Knn <- read.csv("Predict_Knn.csv", header  = TRUE)
row.names(Predict_Knn) <- Predict_Knn$X
Predict_Knn$X <- NULL

Predict_lasso <- read.csv("Predict_lasso.csv", header  = TRUE)
row.names(Predict_lasso) <- Predict_lasso$X
Predict_lasso$X <- NULL

Predict_lda <- read.csv("Predict_lda.csv", header  = TRUE)
row.names(Predict_lda) <- Predict_lda$X
Predict_lda$X <- NULL

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
name_list1 <- my_levels$type...126
name_list2 <- my_levels$C_T...125
name_list3 <- my_levels$Cell_type...128
name_list4 <- my_levels$Cell_type...128[1:4]


#### Predicted XGB ####

Predict_xgb <- Predict_xgb %>%
  mutate(type...126 = case_when(
    type...126 == 1 ~ name_list1[1],
    type...126 == 2 ~ name_list1[2]
  ))


Predict_xgb <- Predict_xgb %>%
  mutate(C_T...125 = case_when(
    C_T...125 == 1 ~ name_list2[1],
    C_T...125 == 2 ~ name_list2[2]
  ))


Predict_xgb <- Predict_xgb %>%
  mutate(Cell_type...128 = case_when(
    Cell_type...128 == 1 ~ name_list3[1],
    Cell_type...128 == 2 ~ name_list3[2],
    Cell_type...128 == 3 ~ name_list3[3],
    Cell_type...128 == 4 ~ name_list3[4],
    Cell_type...128 == 5 ~ name_list3[5],  
  ))


Predict_xgb <- Predict_xgb %>%
  mutate(predicted_xgb = case_when(
    predicted_xgb == 1 ~ name_list4[1],
    predicted_xgb == 2 ~ name_list4[2],
    predicted_xgb == 3 ~ name_list4[3],
    predicted_xgb == 4 ~ name_list4[4] 
  ))
view(Predict_xgb)



#### Predicted KNN ####

Predict_Knn <- Predict_Knn %>%
  mutate(type...126 = case_when(
    type...126 == 1 ~ name_list1[1],
    type...126 == 2 ~ name_list1[2]
  ))


Predict_Knn <- Predict_Knn %>%
  mutate(C_T...125 = case_when(
    C_T...125 == 1 ~ name_list2[1],
    C_T...125 == 2 ~ name_list2[2]
  ))


Predict_Knn <- Predict_Knn %>%
  mutate(Cell_type...128 = case_when(
    Cell_type...128 == 1 ~ name_list3[1],
    Cell_type...128 == 2 ~ name_list3[2],
    Cell_type...128 == 3 ~ name_list3[3],
    Cell_type...128 == 4 ~ name_list3[4],
    Cell_type...128 == 5 ~ name_list3[5],  
  ))


Predict_Knn <- Predict_Knn %>%
  mutate(predicted_Knn = case_when(
    predicted_Knn == 1 ~ name_list4[1],
    predicted_Knn == 2 ~ name_list4[2],
    predicted_Knn == 3 ~ name_list4[3],
    predicted_Knn == 4 ~ name_list4[4] 
  ))
view(Predict_Knn)




#### Predicted LASSO ####

Predict_lasso <- Predict_lasso %>%
  mutate(type...126 = case_when(
    type...126 == 1 ~ name_list1[1],
    type...126 == 2 ~ name_list1[2]
  ))


Predict_lasso <- Predict_lasso %>%
  mutate(C_T...125 = case_when(
    C_T...125 == 1 ~ name_list2[1],
    C_T...125 == 2 ~ name_list2[2]
  ))


Predict_lasso <- Predict_lasso %>%
  mutate(Cell_type...128 = case_when(
    Cell_type...128 == 1 ~ name_list3[1],
    Cell_type...128 == 2 ~ name_list3[2],
    Cell_type...128 == 3 ~ name_list3[3],
    Cell_type...128 == 4 ~ name_list3[4],
    Cell_type...128 == 5 ~ name_list3[5],  
  ))


Predict_lasso <- Predict_lasso %>%
  mutate(Predict_LASSO = case_when(
    Predict_LASSO == 1 ~ name_list4[1],
    Predict_LASSO == 2 ~ name_list4[2],
    Predict_LASSO == 3 ~ name_list4[3],
    Predict_LASSO == 4 ~ name_list4[4] 
  ))

view(Predict_lasso)

#### unify the results ####
# Create a new data frame 'All_results' by selecting only the 'Cell_type...128' column from 'my_data'
All_results <- subset(my_data, select = Cell_type...128)  


# Add predicted values from the Random Forest model to 'predicted_RF' column
All_results$predicted_RF <- Predict_rf$predicted_rf

# Add predicted values from the K-Nearest Neighbors model to 'predicted_KNN' column
All_results$predicted_KNN <- Predict_Knn$predicted_Knn

# Add predicted values from the XGBoost model to 'predicted_XGB' column
All_results$predicted_XGB <- Predict_xgb$predicted_xgb

# Add predicted values from the Naive Bayes model to 'predicted_NB' column
All_results$predicted_NB <- Predict_nb$predicted_nb

# Add predicted values from the Linear Discriminant Analysis model to 'predicted_LDA' column
All_results$predicted_LDA <- Predict_lda$predicted_LDA

# Add predicted values from the LASSO model to 'predicted_LASSO' column
All_results$predicted_LASSO <- Predict_lasso$Predict_LASSO



#### Encode and compare with basic functions ####

columns_to_encode <- c('Cell_type...128','predicted_RF', 'predicted_KNN', 'predicted_XGB', 'predicted_NB', 'predicted_LDA', 'predicted_LASSO')

# Convert specified columns to factor type
All_results <- All_results %>% mutate_at(columns_to_encode, as.factor)
setwd("D:/VarieTHOM/University/QCB/3_SEMESTRE/Data Mining/Laboratory (Blanzieri)/0_PROJECT/Datasets_finals")
write.csv(All_results, "All_results.csv", row.names = TRUE)

summary(All_results[, c('Cell_type...128','predicted_RF', 'predicted_KNN', 'predicted_XGB', 'predicted_NB', 'predicted_LDA', 'predicted_LASSO')])

pairs(All_results[, c('Cell_type...128','predicted_RF', 'predicted_KNN', 'predicted_XGB', 'predicted_NB', 'predicted_LDA', 'predicted_LASSO')])


#### Compare with a logical column ####

# Selecting the relevant columns
prediction_columns <- c('predicted_RF', 'predicted_KNN', 'predicted_XGB', 'predicted_NB')

# Check if values are the same across columns for each row
all_results_same_values <- apply(All_results[, prediction_columns], 1, function(row) all(row == row[1]))

# Create a logical column indicating if all predictions are the same for each row
All_results$all_predictions_same <- all_results_same_values

# create a df with only the models of interest and write it as csv
Compared_RF_KNN_XGB_NB <- subset(All_results, select = c('Cell_type...128','predicted_RF', 'predicted_KNN', 'predicted_XGB', 'predicted_NB','all_predictions_same'))  
setwd("D:/VarieTHOM/University/QCB/3_SEMESTRE/Data Mining/Laboratory (Blanzieri)/0_PROJECT/Datasets_finals")
write.csv(Compared_RF_KNN_XGB_NB, "Compared_RF_KNN_XGB_NB.csv", row.names = TRUE)














