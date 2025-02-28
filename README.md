# HCHS-SOL Survival Machine Learning Models

This repository contains the weighted R code implementations of six survival prediction models developed as part of my thesis work on all-cause mortality prediction in Hispanic/Latino populations. All models incorporate sample weights (using the `bghhsub_s2` variable) to properly account for the complex survey design of the HCHS/SOL data. The repository compares traditional statistical methods with advanced machine learning techniques.

## Project Overview

The goal of this project is to improve mortality prediction by exploring and comparing six different modeling approaches. Each weighted model has been implemented with cross-validation to estimate model performance (using Harrell’s C-index) and to compute individual risk scores (defined as 1 minus the predicted survival probability). In addition, each script produces survival probability plots to visualize the cumulative survival curves over time.

## Repository Structure

- **1_Cox_Backward_Selection.R**  
  Implements a weighted Cox proportional hazards model with backward selection.
  - Uses the `survival`, `Hmisc`, and `ggplot2` packages.
  - Performs 5-fold cross-validation to calculate the C-index and risk scores.
  - Applies sample weights via the `bghhsub_s2` variable.
  - Plots the cumulative survival probability for the final model.

- **2_Weibull_Backward_Selection.R**  
  Implements a weighted parametric Weibull survival model with backward selection.
  - Uses `survival`, `flexsurv`, `MASS`, `Hmisc`, and `ggplot2`.
  - Conducts 5-fold cross-validation and computes predictions at the maximum observed time.
  - Incorporates sample weights to adjust for the complex survey design.
  - Plots the survival curve based on the fitted Weibull model.

- **3_Survival_Tree.R**  
  Implements a weighted nonparametric survival tree.
  - Uses `ranger` (with `num.trees=1`) and `rpart`/`rpart.plot` to visualize the tree structure.
  - Runs 5-fold cross-validation to evaluate model performance and computes risk scores.
  - Applies sample weights during model fitting.
  - Generates survival probability plots for the tree model.

- **4_Random_Survival_Forest.R**  
  Implements a weighted Random Survival Forest model.
  - Utilizes the `ranger` package to build an ensemble of survival trees.
  - Performs 5-fold cross-validation to estimate the C-index and risk scores.
  - Uses sample weights (via `bghhsub_s2`) in the model fitting process.
  - Fits a final model on the entire dataset and plots the cumulative survival curve.

- **5_Lasso_Cox.R**  
  Implements a weighted Lasso-regularized Cox model using penalized regression.
  - Uses `glmnet` along with `survival` and `Hmisc` for model fitting and evaluation.
  - Employs 5-fold cross-validation to select the optimal lambda, compute non-zero coefficients, and assess performance.
  - Incorporates sample weights during model fitting.
  - Produces a final survival curve plot after fitting the model on the full dataset.

- **6_Gradient_Boosting_Cox.R**  
  Implements a weighted Gradient Boosting Machine (GBM) for survival analysis (Cox model).
  - Uses the `gbm` package along with supporting packages (`dplyr`, `Hmisc`, and `ggplot2`).
  - Conducts 5-fold cross-validation, determines the optimal number of trees, and computes risk scores.
  - Incorporates sample weights to adjust predictions for the complex survey design.
  - Produces survival probability curves and additional combined plots comparing models (including plots of C-index and risk score across models).

## Prerequisites

Before running the scripts, ensure you have the following:

- **R** (version 3.6.0 or later recommended)
- Required R packages:  
  `survival`, `Hmisc`, `ggplot2`, `flexsurv`, `MASS`, `ranger`, `dplyr`, `rpart`, `rpart.plot`, `glmnet`, `gbm`

You can install the missing packages using:

```r
install.packages(c("survival", "Hmisc", "ggplot2", "flexsurv", "MASS", "ranger", "dplyr", "rpart", "rpart.plot", "glmnet", "gbm"))
```

## Data

Each script loads the sample dataset from the CSV file `sample.csv`.

The dataset is preprocessed by converting key categorical variables (e.g., `hisp_strat`, `age_strat`, and several predictors) into factors.

## How to Run

Each R script is self-contained. To run any script:

1. Open the desired script (e.g., `1_Cox_Backward_Selection.R`) in RStudio or your preferred R environment.
2. Ensure that the dataset path in the script is correctly specified.
3. Execute the script. Each script will:
   - Perform 5-fold cross-validation,
   - Output evaluation metrics (mean and standard deviation of the C-index and risk scores),
   - Generate survival probability plots.

You can run the scripts one by one to compare model performance and visualizations.

## Evaluation Metrics

- **C-index (Harrell’s Concordance Index):** Measures the model’s ability to correctly rank survival times.
- **Risk Score:** Calculated as 1 minus the predicted survival probability for each individual.
- Cross-validation (5-fold) is used to provide robust performance estimates.

## Visualization

Each script produces survival probability plots using `ggplot2`. Additionally, the GBM script (`6_Gradient_Boosting_Cox.R`) includes code to generate combined plots comparing:

- Survival curves from all models,
- Mean C-index and risk scores across models.

## Customization

- **Model Parameters:** You can adjust parameters (e.g., number of trees, interaction depth, learning rate for GBM; mtry, number of trees for random forests) directly in the script.
- **Cross-Validation:** The number of folds is set to 5 by default; modify the `folds` variable if needed.
- **Plotting:** Customize plot titles, axis labels, and themes by editing the `ggplot2` sections.

## Acknowledgments

This work is part of my thesis research on improving all-cause mortality prediction in Hispanic/Latino populations using advanced machine learning methods. Special thanks to the HCHS/SOL team for providing the dataset and to my thesis advisor Dr Jianwen Cai and PhD candidate Haolin Li for their support.
