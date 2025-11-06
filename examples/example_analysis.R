# Load required libraries
library(glmnet)
library(randomForest)
library(caret)
library(ggplot2)

# Source the functions
source("R/variable_selection.R")
source("R/preprocessing.R")
source("R/utils.R")

# Set random seed for reproducibility
set.seed(123)

cat("=== Gene Expression Variable Selection Example ===\n\n")

# ============================================
# Step 1: Generate Synthetic Gene Expression Data
# ============================================
cat("Step 1: Generating synthetic gene expression data...\n")

n_samples <- 200    # Number of samples
n_genes <- 1000     # Number of genes
n_informative <- 20 # Number of truly informative genes

# Generate random gene expression data
gene_expression <- matrix(rnorm(n_samples * n_genes), nrow = n_samples, ncol = n_genes)
colnames(gene_expression) <- paste0("Gene_", 1:n_genes)
rownames(gene_expression) <- paste0("Sample_", 1:n_samples)

# Create response variable influenced by subset of genes
true_genes <- paste0("Gene_", 1:n_informative)
true_coefs <- rnorm(n_informative, mean = 0, sd = 2)

phenotype <- gene_expression[, 1:n_informative] %*% true_coefs + rnorm(n_samples, sd = 1)
phenotype <- as.numeric(phenotype)

cat(sprintf("Generated data: %d samples, %d genes\n", n_samples, n_genes))
cat(sprintf("True informative genes: %d\n\n", n_informative))

# ============================================
# Step 2: Data Preprocessing
# ============================================
cat("Step 2: Preprocessing data...\n")

# Add some missing values (10% random)
n_missing <- floor(0.10 * n_samples * n_genes)
missing_idx <- sample(1:(n_samples * n_genes), n_missing)
gene_expression[missing_idx] <- NA

# Run preprocessing pipeline
gene_expression_clean <- preprocess_pipeline(
  gene_expression,
  normalize = TRUE,
  norm_method = "zscore",
  filter_variance = TRUE,
  var_threshold = 0.01,
  handle_na = TRUE,
  na_method = "mean",
  remove_outliers_flag = FALSE,
  remove_corr = TRUE,
  cor_threshold = 0.95
)

cat("\n")

# ============================================
# Step 3: Train-Test Split
# ============================================
cat("Step 3: Splitting data into train and test sets...\n")

split_data <- train_test_split(
  gene_expression_clean,
  phenotype,
  train_ratio = 0.8,
  seed = 123
)

cat("\n")

# ============================================
# Step 4: Variable Selection - LASSO
# ============================================
cat("Step 4a: Running LASSO variable selection...\n")

lasso_result <- lasso_selection(
  split_data$x_train,
  split_data$y_train,
  alpha = 1,
  nfolds = 10
)

cat(sprintf("\nLASSO selected %d genes\n", length(lasso_result$selected_vars_1se)))

# Check overlap with true genes
true_genes_found_lasso <- intersect(lasso_result$selected_vars_1se, true_genes)
cat(sprintf("True informative genes found: %d out of %d\n\n",
            length(true_genes_found_lasso), n_informative))

# ============================================
# Step 5: Variable Selection - Elastic Net
# ============================================
cat("Step 4b: Running Elastic Net variable selection...\n")

elasticnet_result <- elastic_net_selection(
  split_data$x_train,
  split_data$y_train,
  alpha = 0.5,
  nfolds = 10
)

cat(sprintf("\nElastic Net selected %d genes\n", length(elasticnet_result$selected_vars_1se)))

true_genes_found_enet <- intersect(elasticnet_result$selected_vars_1se, true_genes)
cat(sprintf("True informative genes found: %d out of %d\n\n",
            length(true_genes_found_enet), n_informative))

# ============================================
# Step 6: Variable Selection - Random Forest
# ============================================
cat("Step 4c: Running Random Forest variable importance...\n")

rf_result <- rf_variable_importance(
  split_data$x_train,
  split_data$y_train,
  ntree = 500,
  top_n = 50
)

cat(sprintf("\nRandom Forest selected %d genes\n", length(rf_result$selected_vars)))

true_genes_found_rf <- intersect(rf_result$selected_vars, true_genes)
cat(sprintf("True informative genes found: %d out of %d\n\n",
            length(true_genes_found_rf), n_informative))

# ============================================
# Step 7: Model Validation
# ============================================
cat("Step 5: Validating selected variables on test set...\n\n")

# Validate LASSO selection
lasso_validation <- validate_selection(
  split_data$x_train,
  split_data$y_train,
  split_data$x_test,
  split_data$y_test,
  lasso_result$selected_vars_1se,
  method = "lm"
)
print_performance(lasso_validation$metrics, "LASSO")

# Validate Elastic Net selection
enet_validation <- validate_selection(
  split_data$x_train,
  split_data$y_train,
  split_data$x_test,
  split_data$y_test,
  elasticnet_result$selected_vars_1se,
  method = "lm"
)
print_performance(enet_validation$metrics, "Elastic Net")

# Validate Random Forest selection
rf_validation <- validate_selection(
  split_data$x_train,
  split_data$y_train,
  split_data$x_test,
  split_data$y_test,
  rf_result$selected_vars,
  method = "lm"
)
print_performance(rf_validation$metrics, "Random Forest")

# ============================================
# Step 8: Visualization
# ============================================
cat("Step 6: Creating visualizations...\n")

# Plot LASSO path
if (require(graphics)) {
  cat("Plotting LASSO regularization path...\n")
  pdf("examples/lasso_path.pdf", width = 12, height = 6)
  plot_lasso_path(lasso_result)
  dev.off()
}

# Plot Random Forest variable importance
if (require(ggplot2)) {
  cat("Plotting Random Forest variable importance...\n")
  p <- plot_variable_importance(rf_result$importance, top_n = 30,
                                title = "Top 30 Genes by Random Forest Importance")
  ggsave("examples/rf_importance.pdf", p, width = 10, height = 8)
}

cat("\n")

# ============================================
# Step 9: Export Results
# ============================================
cat("Step 7: Exporting results...\n")

# Export selected genes
export_selected_vars(
  lasso_result$selected_vars_1se,
  "examples/lasso_selected_genes.csv",
  method_name = "LASSO"
)

export_selected_vars(
  rf_result$selected_vars,
  "examples/rf_selected_genes.csv",
  method_name = "Random Forest"
)

# Generate reports
generate_report(lasso_result, "LASSO", "examples/lasso_report.txt")
generate_report(rf_result, "Random Forest", "examples/rf_report.txt")

# Create gene set enrichment input
create_enrichment_input(
  lasso_result$selected_vars_1se,
  "examples/genes_for_enrichment.txt"
)

cat("\n")

# ============================================
# Step 10: Summary
# ============================================
cat("=== Analysis Summary ===\n")
cat(sprintf("Dataset: %d samples, %d genes (after preprocessing)\n",
            nrow(gene_expression_clean), ncol(gene_expression_clean)))
cat(sprintf("Training set: %d samples\n", nrow(split_data$x_train)))
cat(sprintf("Test set: %d samples\n\n", nrow(split_data$x_test)))

cat("Variables selected by each method:\n")
cat(sprintf("  LASSO:       %d genes (R² = %.3f)\n",
            length(lasso_result$selected_vars_1se), lasso_validation$metrics$R_squared))
cat(sprintf("  Elastic Net: %d genes (R² = %.3f)\n",
            length(elasticnet_result$selected_vars_1se), enet_validation$metrics$R_squared))
cat(sprintf("  Random Forest: %d genes (R² = %.3f)\n\n",
            length(rf_result$selected_vars), rf_validation$metrics$R_squared))

cat("True informative genes recovered:\n")
cat(sprintf("  LASSO:         %d / %d (%.1f%%)\n",
            length(true_genes_found_lasso), n_informative,
            100 * length(true_genes_found_lasso) / n_informative))
cat(sprintf("  Elastic Net:   %d / %d (%.1f%%)\n",
            length(true_genes_found_enet), n_informative,
            100 * length(true_genes_found_enet) / n_informative))
cat(sprintf("  Random Forest: %d / %d (%.1f%%)\n",
            length(true_genes_found_rf), n_informative,
            100 * length(true_genes_found_rf) / n_informative))

cat("\n=== Analysis Complete ===\n")
cat("Results have been saved to the 'examples' directory\n")
