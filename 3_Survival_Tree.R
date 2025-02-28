# Load necessary libraries
library(ranger)
library(survival)
library(Hmisc)
library(dplyr)

# Load the dataset
data <- read.csv("/Users/liuxingyan/Desktop/HCHS:SOL Codes/HCHS survival data/sample.csv")
columns_to_factor <- c("hisp_strat", "age_strat", "x1", "x2", "x3", "x4","x12","x13","x14","x18")
for (col in columns_to_factor) {
  data[[col]] <- as.factor(data[[col]])
}

# Define the formula for the Random Survival Forest model
rf_formula <- as.formula(Surv(time, status) ~ hisp_strat + age_strat + x1 + x2 + x3 + x4 + x8 + x12 + x13 + x14 + x15 + x17 + x18)

# Set seed for reproducibility
set.seed(123)

# Initialize variables for 5-fold cross-validation
folds <- 5
cindex_values_survival_tree <- c()
risk_scores_survival_tree <- list()

# Perform 5-fold cross-validation
for (i in 1:folds) {
  cat("Fold:", i, "\n")
  
  # Create train and test indices
  train_indices <- sample(1:nrow(data), size = 0.8 * nrow(data))
  test_indices <- setdiff(1:nrow(data), train_indices)
  
  train_data <- data[train_indices, ]
  test_data <- data[test_indices, ]
  
  # Fit the Random Survival Forest model
  rf_model <- ranger(
    rf_formula,
    data = train_data,
    mtry = 13, # Total number of covariates
    num.trees = 1,
    splitrule = 'extratrees',
    case.weights = train_data$bghhsub_s2,
    respect.unordered.factors = 'order',
    verbose = TRUE
  )
  
  # Make predictions
  predictions <- predict(rf_model, data = test_data)
  surv_predictions <- as.data.frame(predictions$survival)
  time_points <- predictions$unique.death.times
  
  # Function to get the survival probability at the left boundary of the interval
  left_boundary_survival_prob <- function(surv_probs, time_points, target_time) {
    idx <- max(which(time_points <= target_time))
    surv_probs[idx]
  }
  
  # Define the target time for prediction
  max_time <- max(test_data$time)
  
  # Apply the function to each patient's survival probabilities
  predicted_survival_at_time <- apply(
    surv_predictions, 
    1, 
    left_boundary_survival_prob, 
    time_points = time_points, 
    target_time = max_time
  )
  
  # Calculate the risk score
  risk_score_survival_tree <- 1 - predicted_survival_at_time
  risk_scores_survival_tree <- c(risk_scores_survival_tree, list(risk_score_survival_tree))
  
  # Calculate the C-index on the test set
  cindex_survival_tree <- rcorr.cens(predicted_survival_at_time, Surv(test_data$time, test_data$status))["C Index"]
  cindex_values_survival_tree <- c(cindex_values_survival_tree, cindex_survival_tree)
}

# Compute the mean and standard deviation of the C-index
mean_cindex_survival_tree <- mean(cindex_values_survival_tree)
std_cindex_survival_tree <- sd(cindex_values_survival_tree)

# Compute the mean and standard deviation of the risk score
risk_scores_survival_tree <- unlist(risk_scores_survival_tree)
mean_risk_score_survival_tree <- mean(risk_scores_survival_tree)
std_risk_score_survival_tree <- sd(risk_scores_survival_tree)

# Print the results
cat("Mean C-index for Survival Tree: ", mean_cindex_survival_tree, "\n")
cat("Standard Deviation of C-index for Survival Tree: ", std_cindex_survival_tree, "\n")
cat("Mean Risk Score for Survival Tree: ", mean_risk_score_survival_tree, "\n")
cat("Standard Deviation of Risk Score for Survival Tree: ", std_risk_score_survival_tree, "\n")




# Fit the final Random Survival Forest model on the entire dataset
final_tree_model <- ranger(
  rf_formula,
  data = data,
  mtry = 13, # Total number of covariates
  num.trees = 1,
  splitrule = 'extratrees',
  case.weights = data$bghhsub_s2,
  respect.unordered.factors = 'order',
  verbose = TRUE
)

# Make predictions on the entire dataset
final_predictions <- predict(final_tree_model, data = data)
final_surv_predictions <- as.data.frame(final_predictions$survival)
final_time_points <- final_predictions$unique.death.times

tree_plot_data <- data.frame(time = final_time_points, surv_prob_tree = colMeans(final_surv_predictions))

# Plot using ggplot2
ggplot(tree_plot_data, aes(x = time, y = surv_prob_tree)) +
  geom_step() +
  labs(title = "Cumulative Survival Probability - Survival Tree", 
       x = "Follow-up Time in Days", 
       y = "Survival Probability") +
  theme_minimal()





library(rpart)
library(rpart.plot)

survival_tree <- rpart(Surv(time, status) ~ hisp_strat + age_strat + x1 + x2 + x3 + x4 + x8 + x12 + x13 + x14 + x15 + x17 + x18,
                       data = data, 
                       weights = data$bghhsub_s2)
rpart.plot(survival_tree, type = 3,fallen.leaves = TRUE, 
           box.palette = "auto", shadow.col = "gray", nn = TRUE)
