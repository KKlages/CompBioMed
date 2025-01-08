install.packages("glmnet")

library(glmnet)

otus = readRDS("OTU_data_spike.rds")
otus_df <- as.data.frame(otus)

indoxyl = read.csv("indoxylSulfate.csv")

indoxyl$X == (rownames(otus))
 

model <- lm(indoxyl$IS3_corrected ~ ., data = otus_df)

predictions <- predict(model)
correlation <- cor(predictions, indoxyl$IS3_corrected)
print(correlation)

predictions_loocv <- numeric(length = nrow(otus_df))
actual_values <- indoxyl$IS3_corrected

for(i in 1:nrow(otus_df)){
  train_data = otus_df[-i,]
  train_response = indoxyl$IS3_corrected[-i]
  model = lm(train_response ~ ., data = train_data)
  test_data <- otus_df[i, ]
  predictions_loocv[i] <- predict(model, newdata = test_data)
}
correlation <- cor(predictions_loocv, actual_values)
print(correlation)

X <- as.matrix(otus_df)
y <- indoxyl$IS3_corrected
cv_fit <- cv.glmnet(X, y, alpha = 1)
predictions_lasso <- predict(cv_fit, newx = X, s = "lambda.min")
plot(cv_fit)

predictions_loocv_lasso <- numeric(length = nrow(otus_df))  # Match length to rows in otus_df
actual_values_lasso <- indoxyl$IS3_corrected

for (i in 1:nrow(otus_df)) {
  X <- as.matrix(otus_df[-i, ])  
  Y <- indoxyl$IS3_corrected[-i]  
  cv_fit <- cv.glmnet(X, Y, alpha = 1)  
  
  test_data <- as.matrix(otus_df[i, , drop = FALSE]) 
  predictions_loocv_lasso[i] <- predict(cv_fit, newx = test_data, s = "lambda.min") 
}

correlation_lasso <- cor(predictions_loocv_lasso, actual_values_lasso)
print(correlation_lasso)
#########################################
# Load data
load("mmml_vsn.rda")

# Extract expression matrix and phenotype data
expression_matrix <- exprs(mmml.vsn)
phenotype_data <- pData(mmml.vsn)

# Set random seed for reproducibility
set.seed(123)

# Calculate the number of samples and split for training/testing
n_samples <- ncol(expression_matrix)
n_train <- floor(2/3 * n_samples)
train_indices <- sample(1:n_samples, size = n_train)
test_indices <- setdiff(1:n_samples, train_indices)

# Split expression data
x.training <- expression_matrix[, train_indices]
x.test <- expression_matrix[, test_indices]

# Create binary outcome variables
y.training <- ifelse(phenotype_data[train_indices, "GCBABC"] == "ABC", 1, 0)
y.test <- ifelse(phenotype_data[test_indices, "GCBABC"] == "ABC", 1, 0)

# Verify the data
cat("Training data dimensions:", dim(x.training), "\n")
cat("Length of y.training:", length(y.training), "\n")
cat("Distribution of classes in training:\n")
print(table(y.training))

# Transpose x.training for glmnet
x.training_t <- t(x.training)

# Perform cross-validated logistic regression (Lasso)
library(glmnet)
cv_fit <- cv.glmnet(x = x.training_t, 
                    y = y.training,
                    family = "binomial",
                    alpha = 1,  # Lasso
                    nfolds = 10)

# Print model summary
cat("Optimal lambda:", cv_fit$lambda.min, "\n")
cat("Number of non-zero coefficients:", sum(coef(cv_fit, s = "lambda.min") != 0), "\n")

# Install and load ROCR package
if (!require("ROCR")) install.packages("ROCR")
library(ROCR)

# Get predictions for test data (probabilities)
x.test_t <- t(x.test)
pred_prob <- predict(cv_fit, newx = x.test_t, s = "lambda.min", type = "response")

# Create prediction object for ROCR
pred_obj <- prediction(pred_prob, y.test)

# Calculate performance metrics
roc_perf <- performance(pred_obj, measure = "tpr", x.measure = "fpr")
auc_perf <- performance(pred_obj, measure = "auc")
auc_value <- auc_perf@y.values[[1]]

# Plot ROC curve
plot(roc_perf, 
     main = "ROC Curve for DLBCL Subtype Classification",
     col = "blue", lwd = 2)
abline(a = 0, b = 1, lty = 2, col = "gray")
legend("bottomright", legend = paste("AUC =", round(auc_value, 3)), col = "blue", lwd = 2)

# Calculate and print additional performance metrics
pred_class <- ifelse(pred_prob > 0.5, 1, 0)
accuracy <- mean(pred_class == y.test)
conf_matrix <- table(Predicted = pred_class, Actual = y.test)

cat("Model Performance Metrics:\n")
cat("AUC:", round(auc_value, 3), "\n")
cat("Test Accuracy:", round(accuracy, 3), "\n")
cat("\nConfusion Matrix:\n")
print(conf_matrix)

# Calculate sensitivity and specificity
sensitivity <- conf_matrix[2,2] / sum(conf_matrix[,2])
specificity <- conf_matrix[1,1] / sum(conf_matrix[,1])
cat("Sensitivity:", round(sensitivity, 3), "\n")
cat("Specificity:", round(specificity, 3), "\n")

# Extract and examine coefficients from Lasso model
lasso_coef <- coef(cv_fit, s = "lambda.min")
nonzero_coef <- lasso_coef[lasso_coef != 0]
coef_sorted <- sort(abs(nonzero_coef), decreasing = TRUE)

cat("Top 10 most important features from Lasso:\n")
print(head(coef_sorted, 10))

# Fit Ridge regression (alpha = 0)
cv_fit_ridge <- cv.glmnet(x = x.training_t, 
                          y = y.training,
                          family = "binomial",
                          alpha = 0,  # Ridge
                          nfolds = 10)

# Compare cross-validation plots
par(mfrow = c(1, 2))
plot(cv_fit, main = "Lasso Cross Validation")
plot(cv_fit_ridge, main = "Ridge Cross Validation")
par(mfrow = c(1, 1))

# Get Ridge coefficients
ridge_coef <- coef(cv_fit_ridge, s = "lambda.min")
ridge_coef_sorted <- sort(abs(as.vector(ridge_coef)), decreasing = TRUE)

cat("\nTop 10 most important features from Ridge:\n")
print(head(ridge_coef_sorted, 10))

# Plot regularization paths
lasso_path <- glmnet(x = x.training_t, 
                     y = y.training,
                     family = "binomial",
                     alpha = 1)
ridge_path <- glmnet(x = x.training_t, 
                     y = y.training,
                     family = "binomial",
                     alpha = 0)

par(mfrow = c(1, 2))
plot(lasso_path, xvar = "lambda", main = "Lasso Regularization Path")
plot(ridge_path, xvar = "lambda", main = "Ridge Regularization Path")
par(mfrow = c(1, 1))

# Compare performance on test set (AUC)
ridge_pred_prob <- predict(cv_fit_ridge, newx = x.test_t, s = "lambda.min", type = "response")
lasso_pred_prob <- predict(cv_fit, newx = x.test_t, s = "lambda.min", type = "response")

ridge_pred <- prediction(ridge_pred_prob, y.test)
lasso_pred <- prediction(lasso_pred_prob, y.test)

ridge_auc <- performance(ridge_pred, measure = "auc")@y.values[[1]]
lasso_auc <- performance(lasso_pred, measure = "auc")@y.values[[1]]

cat("\nModel Performance Comparison:\n")
cat("Ridge AUC:", round(ridge_auc, 3), "\n")
cat("Lasso AUC:", round(lasso_auc, 3), "\n")