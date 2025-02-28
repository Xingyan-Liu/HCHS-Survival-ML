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

# Perform 5-fold cross-validation and compute the C-index and risk score
set.seed(123) # For reproducibility
folds <- 5
cindex_values_rf <- c()
risk_scores_rf <- list()

for (i in 1:folds) {
  cat(i, "\n")
  # Create train and test indices
  train_indices <- sample(1:nrow(data), size = 0.8 * nrow(data))
  test_indices <- setdiff(1:nrow(data), train_indices)
  
  train_data <- data[train_indices, ]
  test_data <- data[test_indices, ]
  
  # Fit the Random Survival Forest model
  rf_model <- ranger(
    rf_formula,
    data = train_data,
    #importance = 'permutation',
    splitrule = 'extratrees',
    case.weights = train_data$bghhsub_s2,
    respect.unordered.factors = 'order',
    verbose = TRUE
  )
  
  # Make predictions
  predictions <- predict(rf_model, data = test_data)
  surv_predictions <- as.data.frame(predictions$survival)
  time_points <- predictions$unique.death.times
  
  # Define a function to get the survival probability at the left boundary of the interval
  left_boundary_survival_prob <- function(surv_probs, time_points, target_time) {
    idx <- max(which(time_points <= target_time))
    surv_probs[idx]
  }
  
  # Define the target time for prediction
  max_time <- max(test_data$time)
  
  # Apply the left boundary function to each patient's survival probabilities
  predicted_survival_at_time <- apply(
    surv_predictions, 
    1, 
    left_boundary_survival_prob, 
    time_points = time_points, 
    target_time = max_time
  )
  
  # Calculate the risk score as 1 - predicted survival probability for each individual
  risk_score_rf <- 1 - predicted_survival_at_time
  risk_scores_rf <- c(risk_scores_rf, list(risk_score_rf))
  
  # Calculate the C-index on the test set
  cindex_rf <- rcorr.cens(predicted_survival_at_time, Surv(test_data$time, test_data$status))["C Index"]
  cindex_values_rf <- c(cindex_values_rf, cindex_rf)
}

# Compute the mean and standard deviation of the C-index
mean_cindex_rf <- mean(cindex_values_rf)
std_cindex_rf <- sd(cindex_values_rf)

# Compute the mean and standard deviation of the risk score
risk_scores_rf <- unlist(risk_scores_rf)
mean_risk_score_rf <- mean(risk_scores_rf)
std_risk_score_rf <- sd(risk_scores_rf)

# Print the results
cat("Mean C-index for Random Survival Forest: ", mean_cindex_rf, "\n")
cat("Standard Deviation of C-index for Random Survival Forest: ", std_cindex_rf, "\n")
cat("Mean Risk Score for Random Survival Forest: ", mean_risk_score_rf, "\n")
cat("Standard Deviation of Risk Score for Random Survival Forest: ", std_risk_score_rf, "\n")




# Fit the final Random Survival Forest model on the entire dataset
final_rf_model <- ranger(
  rf_formula,
  data = data,
  splitrule = 'extratrees',
  case.weights = data$bghhsub_s2,
  respect.unordered.factors = 'order',
  verbose = TRUE
)

# Make predictions on the entire dataset
final_predictions <- predict(final_rf_model, data = data)
final_surv_predictions <- as.data.frame(final_predictions$survival)
final_time_points <- final_predictions$unique.death.times

rf_plot_data <- data.frame(time = final_time_points, surv_prob_rf = colMeans(final_surv_predictions))

# Plot using ggplot2
ggplot(rf_plot_data, aes(x = time, y = surv_prob_rf)) +
  geom_step() +
  labs(title = "Cumulative Survival Probability - Random Survival Forest", 
       x = "Follow-up Time in Days", 
       y = "Survival Probability") +
  theme_minimal()


