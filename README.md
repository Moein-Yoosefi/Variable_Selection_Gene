# Variable Selection for Gene Expression Data

A comprehensive R package for performing variable selection on gene expression and genomic data using multiple state-of-the-art algorithms.

## Overview

This project provides implementations of various variable selection methods specifically designed for high-dimensional gene expression data, where the number of genes (features) often exceeds the number of samples. The package includes preprocessing pipelines, multiple selection algorithms, validation tools, and visualization utilities.

## Features

### Variable Selection Methods

- **LASSO Regression** - L1 regularization for sparse feature selection
- **Elastic Net** - Combined L1 and L2 regularization
- **Random Forest Variable Importance** - Tree-based feature ranking
- **Stepwise Selection** - Forward, backward, or bidirectional selection
- **Boruta Algorithm** - All-relevant feature selection (optional)

### Data Preprocessing

- Missing value imputation (mean, median, KNN)
- Variance filtering for low-variance genes
- Normalization methods (z-score, min-max, quantile, log2)
- Outlier detection and removal
- Correlation-based feature filtering
- Complete preprocessing pipeline

### Utilities

- Model validation and cross-validation
- Performance metrics (RMSE, MAE, R², correlation)
- Visualization tools (importance plots, heatmaps, LASSO paths)
- Export functionality for downstream analysis
- Gene set enrichment analysis input generation

## Installation

### Required R Packages

Install the required packages in R:

```r
# Core packages
install.packages(c("glmnet", "randomForest", "caret", "MASS", "ggplot2", "pheatmap"))

# Optional packages
install.packages("Boruta")

# Bioconductor packages (optional, for advanced normalization)
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("preprocessCore")
BiocManager::install("impute")
```

### Clone Repository

```bash
git clone https://github.com/Moein-Yoosefi/Variable_Selection_Gene.git
cd Variable_Selection_Gene
```

## Project Structure

```
Variable_Selection_Gene/
├── R/
│   ├── variable_selection.R  # Main variable selection algorithms
│   ├── preprocessing.R        # Data preprocessing functions
│   └── utils.R               # Utility and visualization functions
├── examples/
│   └── example_analysis.R    # Complete example workflow
├── data/                     # Place your data files here
├── docs/                     # Documentation
├── scripts/                  # Additional analysis scripts
├── README.md                 # This file
└── REQUIREMENTS.txt          # R package dependencies
```

## Quick Start

### 1. Load the Functions

```r
source("R/variable_selection.R")
source("R/preprocessing.R")
source("R/utils.R")
```

### 2. Prepare Your Data

```r
# Load your gene expression data
# Format: rows = samples, columns = genes
gene_data <- read.csv("your_gene_expression.csv", row.names = 1)
phenotype <- read.csv("your_phenotype.csv")$response

# Preprocess the data
gene_data_clean <- preprocess_pipeline(
  gene_data,
  normalize = TRUE,
  norm_method = "zscore",
  filter_variance = TRUE,
  handle_na = TRUE
)
```

### 3. Split Data

```r
# Split into training and testing sets
split <- train_test_split(gene_data_clean, phenotype, train_ratio = 0.8)
```

### 4. Perform Variable Selection

```r
# LASSO selection
lasso_result <- lasso_selection(split$x_train, split$y_train)

# Random Forest importance
rf_result <- rf_variable_importance(split$x_train, split$y_train, top_n = 50)

# Elastic Net
enet_result <- elastic_net_selection(split$x_train, split$y_train, alpha = 0.5)
```

### 5. Validate Results

```r
# Validate on test set
validation <- validate_selection(
  split$x_train, split$y_train,
  split$x_test, split$y_test,
  lasso_result$selected_vars_1se
)

print_performance(validation$metrics, "LASSO")
```

## Detailed Usage

### LASSO Variable Selection

```r
# Perform LASSO with custom parameters
lasso_result <- lasso_selection(
  x = gene_expression_matrix,
  y = phenotype_vector,
  alpha = 1,        # 1 for LASSO, 0 for Ridge
  nfolds = 10       # Number of CV folds
)

# Access results
selected_genes <- lasso_result$selected_vars_1se
coefficients <- lasso_result$coefficients_1se
optimal_lambda <- lasso_result$lambda_1se
```

### Elastic Net Selection

```r
# Elastic Net combines L1 and L2 regularization
enet_result <- elastic_net_selection(
  x = gene_expression_matrix,
  y = phenotype_vector,
  alpha = 0.5,      # 0.5 balances LASSO and Ridge
  nfolds = 10
)
```

### Random Forest Variable Importance

```r
# Rank genes by importance
rf_result <- rf_variable_importance(
  x = gene_expression_matrix,
  y = phenotype_vector,
  ntree = 500,
  top_n = 50        # Select top 50 genes
)

# View importance scores
head(rf_result$importance)

# Plot importance
p <- plot_variable_importance(rf_result$importance, top_n = 30)
print(p)
```

### Stepwise Selection

```r
# Stepwise variable selection
step_result <- stepwise_selection(
  x = gene_expression_matrix,
  y = phenotype_vector,
  direction = "both",  # "forward", "backward", or "both"
  max_vars = 20        # Maximum variables to select
)
```

### Compare Multiple Methods

```r
# Compare all methods
comparison <- compare_methods(
  x = gene_expression_matrix,
  y = phenotype_vector,
  methods = c("lasso", "elasticnet", "rf", "stepwise")
)

# Evaluate performance
performance_df <- compare_performance(
  comparison,
  x_train, y_train,
  x_test, y_test
)
```

## Preprocessing Pipeline

### Complete Pipeline

```r
processed_data <- preprocess_pipeline(
  gene_matrix,
  normalize = TRUE,
  norm_method = "zscore",          # "zscore", "minmax", "quantile", "log2"
  filter_variance = TRUE,
  var_threshold = 0.01,
  handle_na = TRUE,
  na_method = "mean",              # "mean", "median", "knn", "remove"
  remove_outliers_flag = FALSE,
  remove_corr = TRUE,
  cor_threshold = 0.9
)
```

### Individual Preprocessing Steps

```r
# Missing value imputation
data <- handle_missing_values(data, method = "mean")

# Normalization
data <- normalize_genes(data, method = "zscore")

# Filter low variance genes
data <- filter_low_variance(data, var_threshold = 0.01)

# Remove highly correlated features
data <- remove_correlated_features(data, cor_threshold = 0.9)

# Remove outlier samples
outlier_result <- remove_outliers(data, method = "pca", threshold = 3)
data <- outlier_result$data
```

## Visualization

### LASSO Regularization Path

```r
plot_lasso_path(lasso_result)
```

### Variable Importance Plot

```r
p <- plot_variable_importance(rf_result$importance, top_n = 20)
ggsave("importance_plot.pdf", p)
```

### Heatmap of Selected Genes

```r
plot_gene_heatmap(
  gene_expression_matrix,
  selected_genes,
  title = "Expression of Selected Genes"
)
```

## Export Results

### Export Selected Variables

```r
export_selected_vars(
  selected_genes,
  "selected_genes.csv",
  method_name = "LASSO"
)
```

### Generate Report

```r
generate_report(
  lasso_result,
  "LASSO",
  "lasso_analysis_report.txt"
)
```

### Create Gene Set Enrichment Input

```r
# Create input file for tools like DAVID, Enrichr, g:Profiler
create_enrichment_input(
  selected_genes,
  "genes_for_enrichment.txt"
)
```

## Example Analysis

Run the complete example analysis:

```r
source("examples/example_analysis.R")
```

This will:
1. Generate synthetic gene expression data
2. Preprocess the data
3. Split into train/test sets
4. Run multiple variable selection methods
5. Validate results
6. Create visualizations
7. Export results

## Performance Metrics

The package calculates the following metrics:

- **RMSE** - Root Mean Squared Error
- **MAE** - Mean Absolute Error
- **R²** - Coefficient of Determination
- **Correlation** - Pearson correlation between predicted and actual values

## Use Cases

### Cancer Classification

```r
# Select genes that distinguish cancer subtypes
cancer_data <- preprocess_pipeline(tumor_expression_data)
lasso_result <- lasso_selection(cancer_data, cancer_subtype)
```

### Disease Biomarker Discovery

```r
# Identify genes associated with disease status
rf_result <- rf_variable_importance(gene_expression, disease_status, top_n = 100)
biomarker_genes <- rf_result$selected_vars
```

### Survival Analysis

```r
# Select genes predictive of patient survival
enet_result <- elastic_net_selection(gene_data, survival_time, alpha = 0.7)
survival_genes <- enet_result$selected_vars_1se
```

## Tips and Best Practices

1. **Data Preprocessing**: Always preprocess your data before variable selection
2. **Cross-Validation**: Use cross-validation to avoid overfitting
3. **Multiple Methods**: Compare multiple selection methods for robust results
4. **Sample Size**: Ensure adequate sample size (generally n > 30 for reliable results)
5. **Validation**: Always validate selected variables on an independent test set
6. **Biological Validation**: Verify selected genes using biological databases and literature

## Common Issues

### Too Many Variables Selected

- Increase regularization (use lambda.1se instead of lambda.min in LASSO)
- Reduce the number of top variables in Random Forest
- Use stricter variance filtering

### Too Few Variables Selected

- Decrease regularization strength
- Increase top_n parameter
- Lower variance threshold

### Memory Issues with Large Datasets

- Use stepwise selection on a subset of genes
- Filter genes by variance before selection
- Consider using chunked processing

## Citation

If you use this package in your research, please cite:

```
Variable Selection for Gene Expression Data
GitHub: https://github.com/Moein-Yoosefi/Variable_Selection_Gene
```

## References

- Tibshirani, R. (1996). Regression shrinkage and selection via the lasso. Journal of the Royal Statistical Society: Series B (Methodological), 58(1), 267-288.
- Zou, H., & Hastie, T. (2005). Regularization and variable selection via the elastic net. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 67(2), 301-320.
- Breiman, L. (2001). Random forests. Machine learning, 45(1), 5-32.

## License

This project is open source and available for academic and research use.

## Contact

For questions, issues, or contributions, please open an issue on GitHub.

## Contributing

Contributions are welcome! Please feel free to submit pull requests or open issues for bugs and feature requests.
