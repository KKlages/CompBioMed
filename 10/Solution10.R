install.packages("e1071")
library(e1071)
data_pqn <- readRDS("f:/CoBi/BioMed/Tutorial/10/data/data_PQN.rds")
responses <- readRDS("f:/CoBi/BioMed/Tutorial/10/data/responses.rds")

# Create SVM classifier
# Using 'responses' as the target variable and data_pqn as the features
# kernel="radial" is usually a good default choice
# Setting probability=TRUE allows for probability predictions if needed
svm_model_radial <- svm(x = data_pqn, 
                 y = responses, 
                 type = "C-classification", 
                 kernel = "radial",
                 probability = TRUE)

# Print the model summary
print(svm_model)

svm_model_polynomial <- svm(x = data_pqn, 
                 y = responses, 
                 type = "C-classification", 
                 kernel = "polynomial",
                 probability = TRUE)

# Print the model summary
print(svm_model_polynomial)
# Get predictions from both models
predictions_radial <- predict(svm_model_radial, data_pqn, probability = TRUE)
predictions_polynomial <- predict(svm_model_polynomial, data_pqn, probability = TRUE)

# Get probability estimates for each class
prob_radial <- attr(predictions_radial, "probabilities")
prob_polynomial <- attr(predictions_polynomial, "probabilities")

# Print the first few probability predictions
print("Radial Kernel Probabilities:")
head(prob_radial)

print("Polynomial Kernel Probabilities:")
head(prob_polynomial)

# Install and load ROCR package if not already installed
if (!require("ROCR")) {
  install.packages("ROCR")
}
library(ROCR)

# Get probability predictions for both models
predictions_radial <- predict(svm_model_radial, data_pqn, probability = TRUE)
predictions_polynomial <- predict(svm_model_polynomial, data_pqn, probability = TRUE)

# Get probability estimates for each class
prob_radial <- attr(predictions_radial, "probabilities")
prob_polynomial <- attr(predictions_polynomial, "probabilities")

# Create prediction objects for ROCR
pred_radial <- prediction(prob_radial[,2], responses)
pred_polynomial <- prediction(prob_polynomial[,2], responses)

# Calculate performance metrics
perf_radial <- performance(pred_radial, "tpr", "fpr")
perf_polynomial <- performance(pred_polynomial, "tpr", "fpr")

# Calculate AUC
auc_radial <- performance(pred_radial, "auc")@y.values[[1]]
auc_polynomial <- performance(pred_polynomial, "auc")@y.values[[1]]

# Plot ROC curves
plot(perf_radial, col="red", main="ROC Curves for SVM Models")
plot(perf_polynomial, col="blue", add=TRUE)
legend("bottomright", 
       legend=c(paste("Radial (AUC =", round(auc_radial, 3), ")"),
                paste("Polynomial (AUC =", round(auc_polynomial, 3), ")")),
       col=c("red", "blue"), 
       lty=1)

# Function to perform LOOCV with debugging
perform_loocv <- function(data, responses) {
  n <- length(responses)
  predictions_radial <- numeric(n)
  predictions_polynomial <- numeric(n)
  
  # Convert responses to factor if not already
  responses <- as.factor(responses)
  
  # Print data dimensions and summary
  print(paste("Data dimensions:", paste(dim(data), collapse="x")))
  print("Response levels:")
  print(table(responses))
  
  for(i in 1:n) {
    tryCatch({
      # Split data into training and test sets
      train_indices <- setdiff(1:n, i)
      test_index <- i
      
      train_data <- data[train_indices, , drop=FALSE]
      test_data <- data[test_index, , drop=FALSE]
      train_responses <- responses[train_indices]
      
      # Check for NAs in training data
      if(any(is.na(train_data)) || any(is.na(train_responses))) {
        print(paste("Warning: NAs found in training data at iteration", i))
      }
      
      # Train radial model
      svm_model_radial <- svm(x = train_data, 
                             y = train_responses, 
                             type = "C-classification",
                             kernel = "radial",
                             scale = TRUE)  # Added scaling
      
      # Train polynomial model
      svm_model_polynomial <- svm(x = train_data, 
                                 y = train_responses, 
                                 type = "C-classification",
                                 kernel = "polynomial",
                                 scale = TRUE)  # Added scaling
      
      # Make predictions
      pred_radial <- predict(svm_model_radial, test_data)
      pred_polynomial <- predict(svm_model_polynomial, test_data)
      
      # Store predictions
      predictions_radial[i] <- as.character(pred_radial)
      predictions_polynomial[i] <- as.character(pred_polynomial)
      
    }, error = function(e) {
      print(paste("Error at iteration", i, ":", e$message))
      predictions_radial[i] <- NA
      predictions_polynomial[i] <- NA
    })
  }
  
  # Convert predictions to factors
  predictions_radial <- factor(predictions_radial, levels = levels(responses))
  predictions_polynomial <- factor(predictions_polynomial, levels = levels(responses))
  
  # Print summary of predictions
  print("Summary of predictions:")
  print("Radial kernel predictions:")
  print(table(predictions_radial, useNA = "always"))
  print("Polynomial kernel predictions:")
  print(table(predictions_polynomial, useNA = "always"))
  
  return(list(
    radial = predictions_radial,
    polynomial = predictions_polynomial
  ))
}

# Let's first check the data
print("Checking data structure:")
str(data_pqn)
print("Checking responses structure:")
str(responses)

# Remove any rows with NA values
complete_cases <- complete.cases(data_pqn, responses)
data_pqn <- data_pqn[complete_cases, ]
responses <- responses[complete_cases]

# Perform LOOCV
loocv_results <- perform_loocv(data_pqn, responses)

# Calculate accuracy, handling NAs
accuracy_radial <- mean(loocv_results$radial == responses, na.rm = TRUE)
accuracy_polynomial <- mean(loocv_results$polynomial == responses, na.rm = TRUE)

# Print results
print(paste("Radial kernel LOOCV accuracy:", round(accuracy_radial, 3)))
print(paste("Polynomial kernel LOOCV accuracy:", round(accuracy_polynomial, 3)))

# Load required packages
if (!require("ROCR")) {
  install.packages("ROCR")
}
library(ROCR)

# Modify the LOOCV function to return probabilities
perform_loocv_probs <- function(data, responses) {
  n <- length(responses)
  probs_radial <- numeric(n)
  probs_polynomial <- numeric(n)
  
  responses <- as.factor(responses)
  
  for(i in 1:n) {
    train_indices <- setdiff(1:n, i)
    test_index <- i
    
    train_data <- data[train_indices, , drop=FALSE]
    test_data <- data[test_index, , drop=FALSE]
    train_responses <- responses[train_indices]
    
    # Train models
    svm_model_radial <- svm(x = train_data, 
                           y = train_responses, 
                           type = "C-classification",
                           kernel = "radial",
                           probability = TRUE,
                           scale = TRUE)
    
    svm_model_polynomial <- svm(x = train_data, 
                               y = train_responses, 
                               type = "C-classification",
                               kernel = "polynomial",
                               probability = TRUE,
                               scale = TRUE)
    
    # Get probability predictions
    pred_radial <- predict(svm_model_radial, test_data, probability = TRUE)
    pred_polynomial <- predict(svm_model_polynomial, test_data, probability = TRUE)
    
    # Extract probabilities for class 1
    probs_radial[i] <- attr(pred_radial, "probabilities")[, "1"]
    probs_polynomial[i] <- attr(pred_polynomial, "probabilities")[, "1"]
  }
  
  return(list(
    radial = probs_radial,
    polynomial = probs_polynomial
  ))
}

# Get probability predictions
prob_results <- perform_loocv_probs(data_pqn, responses)

# Create prediction objects for ROCR
pred_obj_radial <- prediction(prob_results$radial, responses)
pred_obj_polynomial <- prediction(prob_results$polynomial, responses)

# Calculate ROC curves
roc_radial <- performance(pred_obj_radial, "tpr", "fpr")
roc_polynomial <- performance(pred_obj_polynomial, "tpr", "fpr")

# Calculate PR curves
pr_radial <- performance(pred_obj_radial, "prec", "rec")
pr_polynomial <- performance(pred_obj_polynomial, "prec", "rec")

# Calculate AUCs
auc_roc_radial <- performance(pred_obj_radial, "auc")@y.values[[1]]
auc_roc_polynomial <- performance(pred_obj_polynomial, "auc")@y.values[[1]]

# Calculate AUPRC (Area Under Precision-Recall Curve)
auc_pr_radial <- performance(pred_obj_radial, "aucpr")@y.values[[1]]
auc_pr_polynomial <- performance(pred_obj_polynomial, "aucpr")@y.values[[1]]

# Create ROC plot
par(mfrow=c(1,2))  # Set up 1x2 plotting layout

# Plot ROC curves
plot(roc_radial, 
     col="red", 
     main="ROC Curves",
     lwd=2)
plot(roc_polynomial, 
     col="blue", 
     add=TRUE,
     lwd=2)
abline(a=0, b=1, lty=2, col="gray")  # Add diagonal reference line
legend("bottomright", 
       legend=c(paste("Radial (AUC =", round(auc_roc_radial, 3), ")"),
                paste("Polynomial (AUC =", round(auc_roc_polynomial, 3), ")")),
       col=c("red", "blue"), 
       lwd=2)

# Plot PR curves
plot(pr_radial, 
     col="red", 
     main="Precision-Recall Curves",
     lwd=2)
plot(pr_polynomial, 
     col="blue", 
     add=TRUE,
     lwd=2)
legend("bottomright", 
       legend=c(paste("Radial (AUC =", round(auc_pr_radial, 3), ")"),
                paste("Polynomial (AUC =", round(auc_pr_polynomial, 3), ")")),
       col=c("red", "blue"), 
       lwd=2)

# Print AUC values
cat("\nROC AUC values:\n")
cat("Radial kernel ROC AUC:", round(auc_roc_radial, 3), "\n")
cat("Polynomial kernel ROC AUC:", round(auc_roc_polynomial, 3), "\n")

cat("\nPR AUC values:\n")
cat("Radial kernel PR AUC:", round(auc_pr_radial, 3), "\n")
cat("Polynomial kernel PR AUC:", round(auc_pr_polynomial, 3), "\n")

# Function to perform LOOCV with multiple kernels
perform_multi_kernel_loocv <- function(data, responses) {
  n <- length(responses)
  kernels <- c("radial", "polynomial", "linear", "sigmoid")
  
  # Initialize lists to store results
  probs <- list()
  predictions <- list()
  
  responses <- as.factor(responses)
  
  # Initialize storage for each kernel
  for(kernel in kernels) {
    probs[[kernel]] <- numeric(n)
    predictions[[kernel]] <- numeric(n)
  }
  
  for(i in 1:n) {
    train_indices <- setdiff(1:n, i)
    test_index <- i
    
    train_data <- data[train_indices, , drop=FALSE]
    test_data <- data[test_index, , drop=FALSE]
    train_responses <- responses[train_indices]
    
    # Train and predict with each kernel
    for(kernel in kernels) {
      tryCatch({
        # Train model
        svm_model <- svm(x = train_data, 
                        y = train_responses, 
                        type = "C-classification",
                        kernel = kernel,
                        probability = TRUE,
                        scale = TRUE)
        
        # Get predictions
        pred <- predict(svm_model, test_data, probability = TRUE)
        predictions[[kernel]][i] <- as.numeric(as.character(pred))
        probs[[kernel]][i] <- attr(pred, "probabilities")[, "1"]
      }, error = function(e) {
        cat("Error with kernel", kernel, "at iteration", i, ":", e$message, "\n")
        predictions[[kernel]][i] <- NA
        probs[[kernel]][i] <- NA
      })
    }
  }
  
  return(list(
    probabilities = probs,
    predictions = predictions
  ))
}

# Function to calculate performance metrics
calculate_metrics <- function(predictions, probabilities, true_values) {
  # Convert predictions and true values to factors with same levels
  pred_factor <- factor(predictions, levels = levels(as.factor(true_values)))
  true_factor <- factor(true_values, levels = levels(as.factor(true_values)))
  
  # Create ROCR prediction object
  pred_obj <- prediction(probabilities, true_values)
  
  # Calculate ROC and PR curves
  roc <- performance(pred_obj, "tpr", "fpr")
  pr <- performance(pred_obj, "prec", "rec")
  
  # Calculate AUCs
  auc_roc <- performance(pred_obj, "auc")@y.values[[1]]
  auc_pr <- performance(pred_obj, "aucpr")@y.values[[1]]
  
  # Calculate accuracy
  accuracy <- mean(pred_factor == true_factor, na.rm = TRUE)
  
  # Calculate confusion matrix
  conf_matrix <- table(Predicted = pred_factor, Actual = true_factor)
  
  # Calculate sensitivity (true positive rate)
  sensitivity <- conf_matrix[2,2] / sum(conf_matrix[,2])
  
  # Calculate specificity (true negative rate)
  specificity <- conf_matrix[1,1] / sum(conf_matrix[,1])
  
  return(list(
    accuracy = accuracy,
    sensitivity = sensitivity,
    specificity = specificity,
    auc_roc = auc_roc,
    auc_pr = auc_pr,
    roc = roc,
    pr = pr
  ))
}

# Perform LOOCV for all kernels
results <- perform_multi_kernel_loocv(data_pqn, responses)

# Calculate metrics for each kernel
kernels <- c("radial", "polynomial", "linear", "sigmoid")
metrics <- list()
for(kernel in kernels) {
  metrics[[kernel]] <- calculate_metrics(
    results$predictions[[kernel]], 
    results$probabilities[[kernel]], 
    responses
  )
}

# Create comparison plots
par(mfrow=c(2,2))

# Plot ROC curves
plot(0:1, 0:1, type="n", main="ROC Curves", xlab="False Positive Rate", 
     ylab="True Positive Rate")
colors <- c("red", "blue", "green", "purple")
for(i in seq_along(kernels)) {
  lines(metrics[[kernels[i]]]$roc@x.values[[1]], 
        metrics[[kernels[i]]]$roc@y.values[[1]], 
        col=colors[i], lwd=2)
}
abline(a=0, b=1, lty=2, col="gray")
legend("bottomright", 
       legend=paste(kernels, "(AUC =", 
                   sapply(metrics, function(x) round(x$auc_roc, 3)), ")"),
       col=colors, lwd=2)

# Plot PR curves
plot(0:1, 0:1, type="n", main="Precision-Recall Curves", 
     xlab="Recall", ylab="Precision")
for(i in seq_along(kernels)) {
  lines(metrics[[kernels[i]]]$pr@x.values[[1]], 
        metrics[[kernels[i]]]$pr@y.values[[1]], 
        col=colors[i], lwd=2)
}
legend("bottomright", 
       legend=paste(kernels, "(AUC =", 
                   sapply(metrics, function(x) round(x$auc_pr, 3)), ")"),
       col=colors, lwd=2)

# Create performance comparison barplot
performance_matrix <- sapply(metrics, function(x) 
  c(x$accuracy, x$sensitivity, x$specificity))
rownames(performance_matrix) <- c("Accuracy", "Sensitivity", "Specificity")
barplot(performance_matrix, beside=TRUE, col=colors,
        main="Performance Metrics Comparison",
        legend.text=kernels)

# Print detailed metrics
cat("\nDetailed Performance Metrics:\n")
cat("==========================\n")
for(kernel in kernels) {
  cat("\nKernel:", kernel, "\n")
  cat("Accuracy:", round(metrics[[kernel]]$accuracy, 3), "\n")
  cat("Sensitivity:", round(metrics[[kernel]]$sensitivity, 3), "\n")
  cat("Specificity:", round(metrics[[kernel]]$specificity, 3), "\n")
  cat("ROC AUC:", round(metrics[[kernel]]$auc_roc, 3), "\n")
  cat("PR AUC:", round(metrics[[kernel]]$auc_pr, 3), "\n")
  cat("--------------------------\n")
}

##############################################
install.packages("randomForest")
library(randomForest)
cleveland_data <- readRDS("f:/CoBi/BioMed/Tutorial/10/data/processed_cleveland.rds")
data_pqn <- readRDS("f:/CoBi/BioMed/Tutorial/10/data/data_PQN.rds")
responses <- readRDS("f:/CoBi/BioMed/Tutorial/10/data/responses.rds")

# Check classes of all variables at once
print("Classes of all variables:")
sapply(cleveland_data, class)

# Check individual variables
print("\nChecking a few individual variables:")
class(cleveland_data$age)  # Example for age variable
class(cleveland_data$sex)  # Example for sex variable

# Set seed for reproducibility
set.seed(123)

# Calculate the number of samples for training (75% of total)
n_samples <- nrow(cleveland_data)
n_train <- floor(0.75 * n_samples)

# Create random indices for training set
train_indices <- sample(1:n_samples, n_train)

# Split the data
train_data <- cleveland_data[train_indices, ]
test_data <- cleveland_data[-train_indices, ]

# Print the dimensions of the splits to verify
print("Dimensions of splits:")
print(paste("Training set:", nrow(train_data), "x", ncol(train_data)))
print(paste("Test set:", nrow(test_data), "x", ncol(test_data)))

# Convert the 'num' variable to a factor (if it isn't already)
train_data$num <- as.factor(train_data$num)
test_data$num <- as.factor(test_data$num)

# Train the random forest model
rf_model <- randomForest(num ~ ., 
                        data = train_data,
                        ntree = 9999,  # number of trees
                        importance = TRUE)  # calculate variable importance

# Print the model summary
print(rf_model)
# Create visualization plots
par(mfrow=c(2,2))

# Plot 1: Error rate vs Number of trees
plot(rf_model, 
     main="Error Rate vs Number of Trees",
     lwd=2)

# Plot 2: Variable Importance Plot
varImpPlot(rf_model,
           main="Variable Importance")



# Reset plotting layout
par(mfrow=c(1,1))

# Get feature importance scores
importance_scores <- importance(rf_model)

# Convert to a data frame for easier handling
importance_df <- data.frame(
  Feature = rownames(importance_scores),
  Importance = importance_scores[, "MeanDecreaseGini"]
)

# Sort by importance (descending)
importance_df <- importance_df[order(importance_df$Importance, decreasing = TRUE), ]

# Create barplot
barplot(importance_df$Importance, 
        names.arg = importance_df$Feature,
        main = "Feature Importance in Random Forest",
        xlab = "Features",
        ylab = "Mean Decrease in Gini",
        las = 2,  # Rotate x-axis labels for better readability
        cex.names = 0.7)  # Adjust text size of feature names

# Make predictions on test data
predictions <- predict(rf_model, test_data)

# Create confusion matrix
conf_matrix <- table(Predicted = predictions, Actual = test_data$num)
print("Confusion Matrix:")
print(conf_matrix)

# Calculate accuracy
accuracy <- sum(diag(conf_matrix)) / sum(conf_matrix)
print(paste("Accuracy:", round(accuracy, 3)))

# Calculate sensitivity (true positive rate) for each class
sensitivity <- diag(conf_matrix) / colSums(conf_matrix)
print("\nSensitivity for each class:")
print(round(sensitivity, 3))

# Calculate specificity for each class
specificity <- numeric(nlevels(test_data$num))
for(i in 1:nlevels(test_data$num)) {
    true_neg <- sum(conf_matrix[-i, -i])
    false_pos <- sum(conf_matrix[i, -i])
    specificity[i] <- true_neg / (true_neg + false_pos)
}
names(specificity) <- levels(test_data$num)
print("\nSpecificity for each class:")
print(round(specificity, 3))

# Get probability predictions for each class
prob_predictions <- predict(rf_model, test_data, type = "prob")

# Create plots for each class
par(mfrow=c(2,2))

# Initialize vectors to store AUC values
roc_aucs <- numeric(ncol(prob_predictions))
pr_aucs <- numeric(ncol(prob_predictions))

# Create ROC and PR curves for each class
for(i in 1:ncol(prob_predictions)) {
    # Create binary labels for current class
    binary_labels <- ifelse(test_data$num == levels(test_data$num)[i], 1, 0)
    
    # Create prediction object
    pred <- prediction(prob_predictions[,i], binary_labels)
    
    # ROC curve
    roc_perf <- performance(pred, "tpr", "fpr")
    roc_aucs[i] <- performance(pred, "auc")@y.values[[1]]
    
    # PR curve
    pr_perf <- performance(pred, "prec", "rec")
    pr_aucs[i] <- performance(pred, "aucpr")@y.values[[1]]
    
    # Plot ROC curve
    if(i == 1) {
        plot(roc_perf, 
             main="ROC Curves", 
             col=i,
             lwd=2)
        abline(0, 1, lty=2, col="gray")
    } else {
        plot(roc_perf, 
             add=TRUE, 
             col=i,
             lwd=2)
    }
    
    # Plot PR curve
    if(i == 1) {
        plot(pr_perf, 
             main="Precision-Recall Curves", 
             col=i,
             lwd=2)
    } else {
        plot(pr_perf, 
             add=TRUE, 
             col=i,
             lwd=2)
    }
}

# Add legends
legend("bottomright", 
       legend=paste("Class", levels(test_data$num), 
                   "\nROC AUC =", round(roc_aucs, 3)),
       col=1:ncol(prob_predictions), 
       lwd=2)

# Reset plotting parameters
par(mfrow=c(1,1))

# Print AUC values
cat("\nROC AUC values:\n")
for(i in 1:length(roc_aucs)) {
    cat("Class", levels(test_data$num)[i], ":", round(roc_aucs[i], 3), "\n")
}

cat("\nPR AUC values:\n")
for(i in 1:length(pr_aucs)) {
    cat("Class", levels(test_data$num)[i], ":", round(pr_aucs[i], 3), "\n")
}