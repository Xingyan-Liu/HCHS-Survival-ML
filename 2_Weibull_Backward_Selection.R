# Load necessary libraries
library(survival)
library(MASS)
library(flexsurv)
library(Hmisc)
library(ggplot2)

# Load the dataset
data <- read.csv("/Users/liuxingyan/Desktop/HCHS:SOL Codes/HCHS survival data/sample.csv")
columns_to_factor <- c("hisp_strat", "age_strat", "x1", "x2", "x3", "x4","x12","x13","x14","x18")
for (col in columns_to_factor) {
  data[[col]] <- as.factor(data[[col]])
}

# Define the formula for the Weibull model
weibull_formula <- as.formula(Surv(time, status) ~ hisp_strat + age_strat + x1 + x2 + x3 + x4 + x8 + x12 + x13 + x14 + x15 + x17 + x18)

# Perform 5-fold cross-validation and compute the C-index and risk score
set.seed(123) # For reproducibility
folds <- 5
cindex_values_weibull <- c()
risk_scores_weibull <- list()

predict_survival_weibull <- function(object, newdata, t){
  mu_hat <- predict(object, newdata = newdata, type = "link")
  cum_hazard <- (t / exp(mu_hat))^(1 / object$scale)
  surv <- exp(-cum_hazard)
  return(surv)
}

for (i in 1:folds) {
  cat(i, "\n")
  # Create train and test indices
  train_indices <- sample(1:nrow(data), size = 0.8 * nrow(data))
  test_indices <- setdiff(1:nrow(data), train_indices)
  
  train_data <- data[train_indices, ]
  test_data <- data[test_indices, ]
  
  # Fit the full Weibull model
  full_model <- survreg(weibull_formula, data = train_data, weights = bghhsub_s2, dist = "weibull")
  
  # Perform backward selection using AIC
  backward_model <- step(full_model, direction = "backward", trace = 0)
  
  # Refit the final model after backward selection
  final_model <- survreg(backward_model[["call"]][["formula"]], data = train_data, weights = bghhsub_s2, dist = "weibull")
  
  # Predict the survival probability at the maximum observed time for the test set
  max_time <- max(test_data$time)
  predicted_survival <- predict_survival_weibull(final_model, test_data, t = max_time)
  
  # Calculate the risk score as 1 - predicted survival probability
  risk_score_weibull <- 1 - predicted_survival
  risk_scores_weibull <- c(risk_scores_weibull, list(risk_score_weibull))
  
  # Calculate the C-index on the test set
  cindex_weibull <- rcorr.cens(predicted_survival, Surv(test_data$time, test_data$status))["C Index"]
  cindex_values_weibull <- c(cindex_values_weibull, cindex_weibull)
}

# Compute the mean and standard deviation of the C-index
mean_cindex_weibull <- mean(cindex_values_weibull)
std_cindex_weibull <- sd(cindex_values_weibull)

# Compute the mean and standard deviation of the risk score
risk_scores_weibull <- unlist(risk_scores_weibull)
mean_risk_score_weibull <- mean(risk_scores_weibull)
std_risk_score_weibull <- sd(risk_scores_weibull)

# Print the results
cat("Mean C-index for Weibull: ", mean_cindex_weibull, "\n")
cat("Standard Deviation of C-index for Weibull: ", std_cindex_weibull, "\n")
cat("Mean Risk Score for Weibull: ", mean_risk_score_weibull, "\n")
cat("Standard Deviation of Risk Score for Weibull: ", std_risk_score_weibull, "\n")



# Fit the final Weibull model on the entire dataset
final_weibull_model <- survreg(weibull_formula, data = data, weights = bghhsub_s2, dist = "weibull")

# Define time points for survival probabilities
time_points_weibull <- sort(unique(data$time[data$status == 1]))

# Predict survival probabilities for the entire dataset at each time point
predict_survival_curve <- function(object, newdata, time_points){
  surv_probs <- sapply(time_points, function(t){
    predict_survival_weibull(object, newdata, t)
  })
  return(surv_probs)
}

surv_probs_weibull <- as.data.frame(predict_survival_curve(final_weibull_model, data, time_points_weibull))

# Create a data frame for plotting
weibull_plot_data <- data.frame(time = time_points_weibull, surv_prob_weibull = colMeans(surv_probs_weibull))

# Plot using ggplot2
ggplot(weibull_plot_data, aes(x = time, y = surv_prob_weibull)) +
  geom_step() +
  labs(title = "Cumulative Survival Probability - Weibull Backward Selection", 
       x = "Follow-up Time in Days", 
       y = "Survival Probability") +
  theme_minimal()
