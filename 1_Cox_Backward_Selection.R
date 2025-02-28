# Load necessary libraries
library(survival)
library(Hmisc)
library(ggplot2)

# Load the dataset
data <- read.csv("/Users/liuxingyan/Desktop/HCHS:SOL Codes/HCHS survival data/sample.csv")
columns_to_factor <- c("hisp_strat", "age_strat", "x1", "x2", "x3", "x4","x12","x13","x14","x18")
for (col in columns_to_factor) {
  data[[col]] <- as.factor(data[[col]])
}

# Define the formula for the Cox model
cox_formula <- as.formula(Surv(time, status) ~ hisp_strat + age_strat + x1 + x2 + x3 + x4 + x8 + x12 + x13 + x14 + x15 + x17 + x18)

# Perform 5-fold cross-validation and compute the C-index and risk score
set.seed(123) # For reproducibility
folds <- 5
cindex_values_cox <- c()
risk_scores_cox <- list()

predict_survival_cox <- function(cox_model, newdata, t){
  # Predict the linear predictor for the new data
  linear_predictor <- predict(cox_model, newdata, type = "lp")
  # Calculate the baseline survival function
  base_surv <- survfit(cox_model, newdata = newdata)
  # Find the baseline survival probability at the specified time
  base_times <- summary(base_surv)$time
  base_surv_at_t <- summary(base_surv, times = t)$surv
  # Calculate the survival probability for each observation in the new data
  surv_prob <- base_surv_at_t^exp(linear_predictor)
  return(surv_prob)
}

for (i in 1:folds) {
  cat("Fold:", i, "\n")
  # Create train and test indices
  train_indices <- sample(1:nrow(data), size = 0.8 * nrow(data))
  test_indices <- setdiff(1:nrow(data), train_indices)
  
  train_data <- data[train_indices, ]
  test_data <- data[test_indices, ]
  
  # Fit the full Cox model
  full_model <- coxph(cox_formula, data = train_data, weights = train_data$bghhsub_s2)
  
  # Perform backward selection using AIC
  backward_model <- step(full_model, direction = "backward", trace = 0)
  
  # Refit the final model after backward selection
  final_model <- coxph(backward_model[["call"]][["formula"]], data = train_data, weights = train_data$bghhsub_s2)
  
  # Predict the survival probabilities for the test data
  surv_fit <- survfit(final_model, newdata = test_data)
  surv_probs <- surv_fit[["surv"]][nrow(surv_fit[["surv"]]), ]
  
  # Calculate the risk score as 1 - predicted survival probability for each individual
  risk_score_cox <- 1 - surv_probs
  risk_scores_cox <- c(risk_scores_cox, list(risk_score_cox))
  
  # Calculate the C-index on the test set
  cindex_cox <- rcorr.cens(surv_probs, Surv(test_data$time, test_data$status))["C Index"]
  cindex_values_cox <- c(cindex_values_cox, cindex_cox)
}

# Compute the mean and standard deviation of the C-index
mean_cindex_cox <- mean(cindex_values_cox)
std_cindex_cox <- sd(cindex_values_cox)

# Compute the mean and standard deviation of the risk score
risk_scores_cox <- unlist(risk_scores_cox)
mean_risk_score_cox <- mean(risk_scores_cox)
std_risk_score_cox <- sd(risk_scores_cox)

# Print the results
cat("Mean C-index for Cox Model: ", mean_cindex_cox, "\n")
cat("Standard Deviation of C-index for Cox Model: ", std_cindex_cox, "\n")
cat("Mean Risk Score for Cox Model: ", mean_risk_score_cox, "\n")
cat("Standard Deviation of Risk Score for Cox Model: ", std_risk_score_cox, "\n")



# Fit the final Cox model on the entire dataset
full_cox_model <- coxph(cox_formula, data = data, weights = data$bghhsub_s2)

# Perform backward selection using AIC
backward_cox_model <- step(full_cox_model, direction = "backward", trace = 0)

# Refit the final model after backward selection
final_cox_backward_model <- coxph(backward_cox_model[["call"]][["formula"]], data = data, weights = data$bghhsub_s2)

# Compute survival probabilities for the entire dataset
surv_fit_full <- survfit(final_cox_backward_model, newdata = data)
surv_probs_full <- surv_fit_full[["surv"]]
surv_probs_full <- data.frame(surv_probs_full)

# Create a data frame for plotting
time.interest_all <- sort(unique(data$time[data$status == 1]))
cox_plot_data <- data.frame(time = surv_fit_full[["time"]], surv_prob_cox = rowMeans(surv_probs_full))
cox_plot_data <- cox_plot_data[cox_plot_data$time %in% time.interest_all, ]

# Plot using ggplot2
ggplot(cox_plot_data, aes(x = time, y = surv_prob_cox)) +
  geom_step() +
  labs(title = "Cumulative Survival Probability - Cox Model", 
       x = "Follow-up Time in Days", 
       y = "Survival Probability") +
  theme_minimal()
