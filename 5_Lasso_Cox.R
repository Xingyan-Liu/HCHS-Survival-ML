# Load necessary libraries
library(survival)
library(glmnet)
library(Hmisc)
library(ggplot2)

# Load the dataset
data <- read.csv("/Users/liuxingyan/Desktop/HCHS:SOL Codes/HCHS survival data/sample.csv")

# Factor the specified columns
columns_to_factor <- c("hisp_strat", "age_strat", "x1", "x2", "x3", "x4", "x12", "x13", "x14", "x18")
for (col in columns_to_factor) {
  data[[col]] <- as.factor(data[[col]])
}

# Define the formula for the Cox model (also used for Lasso)
cox_formula <- as.formula(Surv(time, status) ~ hisp_strat + age_strat + x1 + x2 + x3 + x4 + x8 + x12 + x13 + x14 + x15 + x17 + x18)

# Perform 5-fold cross-validation and compute the C-index and risk score
set.seed(123) # For reproducibility
folds <- 5
cindex_values_lasso <- c()
risk_scores_lasso <- list()

for (i in 1:folds) {
  cat("Fold:", i, "\n")
  
  # Create train and test indices
  train_indices <- sample(1:nrow(data), size = 0.8 * nrow(data))
  test_indices <- setdiff(1:nrow(data), train_indices)
  
  train_data <- data[train_indices, ]
  test_data <- data[test_indices, ]
  
  # Define the survival object for the training data
  surv_obj <- Surv(time = train_data$time, event = train_data$status)
  
  # Create the design matrix with one-hot encoding for factors
  x_train <- model.matrix(cox_formula, data = train_data)[, -1] # Remove intercept
  x_test <- model.matrix(cox_formula, data = test_data)[, -1] # Remove intercept
  
  # Fit a Cox Lasso model with cross-validation
  lasso_model <- cv.glmnet(
    x = x_train, 
    y = surv_obj, 
    weights = train_data$bghhsub_s2, 
    standardize = TRUE, 
    family = "cox", 
    alpha = 1
  )
  
  best_lambda <- lasso_model$lambda.min
  
  final_lasso_model <- glmnet(
    x = x_train, 
    y = surv_obj, 
    weights = train_data$bghhsub_s2, 
    standardize = TRUE, 
    family = "cox", 
    alpha = 1, 
    lambda = best_lambda
  )
  
  lasso_coefs <- coef(final_lasso_model)
  non_zero_coefs <- rownames(lasso_coefs)[which(lasso_coefs != 0)]
  
  # Extract non-zero coefficient columns from x_train and x_test
  final_x_train <- x_train[, non_zero_coefs, drop = FALSE]
  final_x_test <- x_test[, non_zero_coefs, drop = FALSE]
  
  # Merge with original train_data and test_data to include necessary columns
  final_train_data <- cbind(train_data[, c("time", "status", "bghhsub_s2")], final_x_train)
  final_test_data <- cbind(test_data[, c("time", "status", "bghhsub_s2")], final_x_test)
  
  # Fit the final Cox model on the training data
  final_cox_model <- coxph(Surv(time, status) ~ ., data = final_train_data, weights = final_train_data$bghhsub_s2)
  
  # Predict the survival probabilities for the test data
  surv_fit <- survfit(final_cox_model, newdata = final_test_data)
  surv_probs <- surv_fit[["surv"]][nrow(surv_fit[["surv"]]), ]
  
  # Calculate the risk score as 1 - survival probability
  risk_score_lasso <- 1 - surv_probs
  risk_scores_lasso <- c(risk_scores_lasso, list(risk_score_lasso))
  
  # Calculate the C-index on the test set
  cindex_lasso <- rcorr.cens(surv_probs, Surv(test_data$time, test_data$status))["C Index"]
  cindex_values_lasso <- c(cindex_values_lasso, cindex_lasso)
}

# Compute the mean and standard deviation of the C-index
mean_cindex_lasso <- mean(cindex_values_lasso)
std_cindex_lasso <- sd(cindex_values_lasso)

# Compute the mean and standard deviation of the risk score
risk_scores_lasso <- unlist(risk_scores_lasso)
mean_risk_score_lasso <- mean(risk_scores_lasso)
std_risk_score_lasso <- sd(risk_scores_lasso)

# Print the results
cat("Mean C-index for Lasso: ", mean_cindex_lasso, "\n")
cat("Standard Deviation of C-index for Lasso: ", std_cindex_lasso, "\n")
cat("Mean Risk Score for Lasso: ", mean_risk_score_lasso, "\n")
cat("Standard Deviation of Risk Score for Lasso: ", std_risk_score_lasso, "\n")



# Fit the final model on the entire dataset
full_x <- model.matrix(cox_formula, data = data)[, -1] # Remove intercept
surv_obj_full <- Surv(time = data$time, event = data$status)

lasso_model <- cv.glmnet(
  x = full_x, 
  y = surv_obj_full, 
  weights = data$bghhsub_s2, 
  standardize = TRUE, 
  family = "cox", 
  alpha = 1
)

best_lambda <- lasso_model$lambda.min

final_lasso_model_full <- glmnet(
  x = full_x, 
  y = surv_obj_full, 
  weights = data$bghhsub_s2, 
  standardize = TRUE, 
  family = "cox", 
  alpha = 1, 
  lambda = best_lambda
)

lasso_coefs_full <- coef(final_lasso_model_full)
non_zero_coefs_full <- rownames(lasso_coefs_full)[which(lasso_coefs_full != 0)]

# Extract non-zero coefficient columns from full_x
final_full_x <- full_x[, non_zero_coefs_full, drop = FALSE]

# Merge with original data to include necessary columns
final_full_data <- cbind(data[, c("time", "status", "bghhsub_s2")], final_full_x)

# Fit the final Cox model on the full data
final_cox_model_full <- coxph(Surv(time, status) ~ ., data = final_full_data, weights = final_full_data$bghhsub_s2)

# Define time points for survival probabilities
time.interest_all <- sort(unique(data$time[data$status == 1]))

# Predict the survival probabilities for the test data
surv_fit <- survfit(final_cox_model_full, newdata = final_full_data)
surv_probs <- surv_fit[["surv"]]

# Create a data frame for plotting
lasso_plot_data <- data.frame(time = surv_fit[["time"]], surv_prob_lasso = rowMeans(surv_probs))
lasso_plot_data <- lasso_plot_data[lasso_plot_data$time %in% time.interest_all, ]

# Plot using ggplot2
ggplot(lasso_plot_data, aes(x = time, y = surv_prob_lasso)) +
  geom_step() +
  labs(title = "Cumulative Survival Probability - Lasso Cox Model", 
       x = "Follow-up Time in Days", 
       y = "Survival Probability") +
  theme_minimal()
