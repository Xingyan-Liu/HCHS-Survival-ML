# Load necessary libraries
library(survival)
library(gbm)
library(dplyr)
library(Hmisc)

# Load the data
data <- read.csv("/Users/liuxingyan/Desktop/HCHS:SOL Codes/HCHS survival data/sample.csv")
columns_to_factor <- c("hisp_strat", "age_strat", "x1", "x2", "x3", "x4","x12","x13","x14","x18")
for (col in columns_to_factor) {
  data[[col]] <- as.factor(data[[col]])
}

# Perform 5-fold cross-validation and compute the C-index and risk score
set.seed(123) # For reproducibility
folds <- 5
cindex_values_gbm <- c()
risk_scores_gbm <- list()

for (i in 1:folds) {
  cat("Fold:", i, "\n")
  
  # Create train and test indices
  train_indices <- sample(1:nrow(data), size = 0.8 * nrow(data))
  test_indices <- setdiff(1:nrow(data), train_indices)
  
  train_data <- data[train_indices, ]
  test_data <- data[test_indices, ]
  
  # Fit the GBM Cox model
  gbm_cox_model <- gbm(Surv(time, status) ~ hisp_strat + age_strat + x1 + x2 + x3 + x4 + x8 + x12 + x13 + x14 + x15 + x17 + x18,
                       data = train_data,
                       distribution = "coxph",
                       weights = bghhsub_s2,
                       n.trees = 100,          # Number of trees
                       interaction.depth = 3,  # Depth of each tree
                       shrinkage = 0.01,       # Learning rate
                       cv.folds = 5,           # Cross-validation folds
                       keep.data = TRUE)       # Keep a copy of the dataset
  
  # Determine the best number of trees
  best_trees <- gbm.perf(gbm_cox_model, method = "cv")
  
  # Predict the log hazard ratios for the training data
  predicted_risk_train <- predict(gbm_cox_model, 
                                  newdata = train_data, 
                                  n.trees = best_trees, 
                                  type = "link")
  
  # Predict the log hazard ratios for the test data
  predicted_risk_test <- predict(gbm_cox_model, 
                                 newdata = test_data, 
                                 n.trees = best_trees)
  
  # Compute the baseline hazard function
  time.interest <- sort(unique(train_data$time[train_data$status==1]))
  baseline_hazard <- basehaz.gbm(t = train_data$time, 
                                 delta = train_data$status, 
                                 f.x = predicted_risk_train, 
                                 t.eval = time.interest, 
                                 cumulative = TRUE)
  
  # Function to get survival probability at the left boundary of the interval
  left_boundary_survival_prob <- function(predicted_risk, baseline_hazard, time.interest, specif.time) {
    idx <- max(which(time.interest <= specif.time))
    if (length(idx) == 0) {
      return(NA) # If specif.time is smaller than any value in time.interest
    }
    exp(-exp(predicted_risk) * baseline_hazard[idx])
  }
  
  # Define the target time for prediction
  max_time <- max(test_data$time)
  
  # Apply the left boundary function to each patient's survival probabilities
  predicted_survival <- sapply(predicted_risk_test, left_boundary_survival_prob, 
                               baseline_hazard = baseline_hazard, 
                               time.interest = time.interest, 
                               specif.time = max_time)
  
  # Calculate the risk score
  risk_score_gbm <- -predicted_risk_test
  risk_scores_gbm <- c(risk_scores_gbm, list(risk_score_gbm))
  
  # Calculate the C-index on the test set
  cindex_gbm <- rcorr.cens(predicted_survival, Surv(test_data$time, test_data$status))["C Index"]
  cindex_values_gbm <- c(cindex_values_gbm, cindex_gbm)
}

# Compute the mean and standard deviation of the C-index
mean_cindex_gbm <- mean(cindex_values_gbm)
std_cindex_gbm <- sd(cindex_values_gbm)

# Compute the mean and standard deviation of the risk score
risk_scores_gbm <- unlist(risk_scores_gbm)
mean_risk_score_gbm <- mean(risk_scores_gbm)
std_risk_score_gbm <- sd(risk_scores_gbm)

# Print the results
cat("Mean C-index for GBM Cox: ", mean_cindex_gbm, "\n")
cat("Standard Deviation of C-index for GBM Cox: ", std_cindex_gbm, "\n")
cat("Mean Risk Score for GBM Cox: ", mean_risk_score_gbm, "\n")
cat("Standard Deviation of Risk Score for GBM Cox: ", std_risk_score_gbm, "\n")




# Fit the final GBM Cox model on the entire dataset
final_gbm_cox_model <- gbm(Surv(time, status) ~ hisp_strat + age_strat + x1 + x2 + x3 + x4 + x8 + x12 + x13 + x14 + x15 + x17 + x18,
                           data = data,
                           distribution = "coxph",
                           weights = bghhsub_s2,
                           n.trees = 100,          
                           interaction.depth = 3,  
                           shrinkage = 0.01,       
                           cv.folds = 5,           
                           keep.data = TRUE)       

best_trees_final <- gbm.perf(final_gbm_cox_model, method = "cv")

# Predict the log hazard ratios for the entire data
predicted_risk_all <- predict(final_gbm_cox_model, 
                              newdata = data, 
                              n.trees = best_trees_final, 
                              type = "link")

# Compute the baseline hazard function
time.interest_all <- sort(unique(data$time[data$status == 1]))
baseline_hazard_all <- basehaz.gbm(t = data$time, 
                                   delta = data$status, 
                                   f.x = predicted_risk_all, 
                                   t.eval = time.interest_all, 
                                   cumulative = TRUE)

# Create a function to get survival probability for all time points
compute_survival_prob <- function(predicted_risk, baseline_hazard, time.interest) {
  surv_probs <- sapply(predicted_risk, function(risk) {
    exp(-exp(risk) * baseline_hazard)
  })
  return(surv_probs)
}

# Compute survival probabilities for the entire dataset for each unique time point
surv_prob_matrix <- compute_survival_prob(predicted_risk_all, baseline_hazard_all, time.interest_all)

# Aggregate the survival probabilities by taking the mean at each time point
mean_surv_prob <- rowMeans(surv_prob_matrix)

# Create a data frame for plotting
gbm_plot_data <- data.frame(time = time.interest_all, surv_prob_gbm = mean_surv_prob)

# Plot using ggplot2
ggplot(gbm_plot_data, aes(x = time, y = surv_prob_gbm)) +
  geom_step() +
  labs(title = "Cumulative Survival Probability", 
       x = "Follow-up Time in Days", 
       y = "1 - Predicted Probability of Incident Event") +
  theme_minimal()









library(ggplot2)
library(gridExtra)
library(cowplot)

# Merge the Lasso and GBM plot data by time
combined_plot_data <- merge(lasso_plot_data, gbm_plot_data, by = "time", all.y = TRUE)
combined_plot_data <- merge(rf_plot_data, combined_plot_data, by = "time", all.y = TRUE)
combined_plot_data <- merge(tree_plot_data, combined_plot_data, by = "time", all.y = TRUE)
combined_plot_data <- merge(weibull_plot_data, combined_plot_data, by = "time", all.y = TRUE)
combined_plot_data <- merge(cox_plot_data, combined_plot_data, by = "time", all.y = TRUE)

# Subset the data where status is 1
status1_data <- data[data$status == 1, c("time", "hisp_strat", "age_strat")]

# Merge hisp_strat and age_strat into combined_plot_data by time
#combined_plot_data <- merge(combined_plot_data, status1_data, by = "time", all.x = TRUE)














combined_plot_data <- read.csv("/Users/liuxingyan/Desktop/HCHS:SOL Codes/Result/combined_plot_data.csv")

# Define the unique combinations of hisp_strat and age_strat
combinations <- unique(status1_data[, c("hisp_strat", "age_strat")])


library(ggplot2)
library(gridExtra)
library(cowplot)

# Function to plot survival data and annotate with max time and survival probability
plot_survival <- function(subset_data, hisp_strat, age_strat) {
  # Convert factors to logical
  hisp_strat <- as.logical(hisp_strat)
  age_strat <- as.logical(age_strat)
  
  # Determine the title based on the conditions
  if (hisp_strat & age_strat) {
    title <- "All-cause mortality for Hispanic more than 45 years old"
  } else if (!hisp_strat & !age_strat) {
    title <- "All-cause mortality for non-Hispanic less than 45 years old"
  } else if (hisp_strat & !age_strat) {
    title <- "All-cause mortality for Hispanic less than 45 years old"
  } else {
    title <- "All-cause mortality for non-Hispanic more than 45 years old"
  }
  
  ggplot(subset_data, aes(x = time)) +
    geom_step(aes(y = surv_prob_cox, color = "Cox Backward Selection")) +
    geom_step(aes(y = surv_prob_weibull, color = "Weibull Backward Selection")) +
    geom_step(aes(y = surv_prob_tree, color = "Survival Tree")) +
    geom_step(aes(y = surv_prob_rf, color = "Random Survival Forest")) +
    geom_step(aes(y = surv_prob_lasso, color = "Lasso Cox")) +
    geom_step(aes(y = surv_prob_gbm, color = "Gradient Boosting")) +
    labs(title = title, 
         x = "Follow-up Time in Years", 
         y = "Predicted Survival Probability") +
    xlim(4, 7.5) + 
    ylim(0.2, 1) + 
    theme_minimal() +
    theme(
      legend.position = "none", # Remove legend from individual plots
      axis.title.y = element_text(size = 10) # Adjust y-axis label size
    )
}

# Prepare a list to store plots
plots <- list()

# Iterate through each combination and create plots
for (i in 1:nrow(combinations)) {
  hisp_strat_val <- as.logical(combinations$hisp_strat[i])
  age_strat_val <- as.logical(combinations$age_strat[i])
  
  # Subset the combined_plot_data by hisp_strat and age_strat
  subset_data <- combined_plot_data[combined_plot_data$hisp_strat == hisp_strat_val &
                                      combined_plot_data$age_strat == age_strat_val, ]
  
  # Create the plot for the current combination and add to the list
  p <- plot_survival(subset_data, age_strat_val, hisp_strat_val)
  
  plots[[i]] <- p
}

# Extract the legend from one plot
legend_plot <- ggplot(combined_plot_data, aes(x = time)) +
  geom_step(aes(y = surv_prob_cox, color = "Cox Backward Selection")) +
  geom_step(aes(y = surv_prob_weibull, color = "Weibull Backward Selection")) +
  geom_step(aes(y = surv_prob_tree, color = "Survival Tree")) +
  geom_step(aes(y = surv_prob_rf, color = "Random Survival Forest")) +
  geom_step(aes(y = surv_prob_lasso, color = "Lasso Cox")) +
  geom_step(aes(y = surv_prob_gbm, color = "Gradient Boosting")) +
  xlim(4, 7.5) + 
  ylim(0.2, 1) + 
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.title.y = element_text(size = 10) # Adjust y-axis label size for legend plot
  )

legend <- get_legend(legend_plot)

# Arrange the plots and legend
plot_grid(
  plot_grid(plotlist = plots, ncol = 2),
  legend,
  ncol = 1,
  rel_heights = c(1, 0.1)
)



# Create the plot for Comparison of C-Index and Risk Score by Model
ggplot(data, aes(x = Model)) + 
  geom_bar(aes(y = Mean_CIndex), stat = "identity", position = "dodge",  
           fill = rgb(173/255, 200/255, 223/255), alpha = 0.7) +
  geom_line(aes(y = Mean_Risk_Score * max(Mean_CIndex) / max(Mean_Risk_Score), group = 1), 
            color = "red", size = 1) +
  geom_point(aes(y = Mean_Risk_Score * max(Mean_CIndex) / max(Mean_Risk_Score)), 
             color = "red", size = 3) +
  scale_y_continuous(
    name = "Mean C-Index",
    sec.axis = sec_axis(~ . * max(data$Mean_Risk_Score) / max(data$Mean_CIndex), name = "Mean Risk Score")
  ) +
  labs(title = "Comparison of C-Index and Risk Score by Model", x = "Model") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

