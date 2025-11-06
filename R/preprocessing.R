# Data Preprocessing Pipeline for Gene Expression Data
# Author: Gene Analysis Team
# Description: Functions for preprocessing and quality control of genomic data

library(preprocessCore)  # For normalization (install via Bioconductor if needed)
library(stats)

#' Load Gene Expression Data
#'
#' Loads gene expression data from CSV or tab-delimited file
#'
#' @param file_path Path to the data file
#' @param sep Separator (default: ",")
#' @param header Whether file has header
#' @return Data frame with gene expression data
load_gene_data <- function(file_path, sep = ",", header = TRUE) {
  cat(sprintf("Loading data from: %s\n", file_path))

  data <- read.table(file_path, sep = sep, header = header, stringsAsFactors = FALSE)

  cat(sprintf("Loaded %d samples and %d genes\n", nrow(data), ncol(data)))

  return(data)
}


#' Remove Low Variance Genes
#'
#' Filters out genes with variance below a threshold
#'
#' @param gene_matrix Matrix of gene expression (samples x genes)
#' @param var_threshold Minimum variance threshold
#' @return Filtered gene matrix
filter_low_variance <- function(gene_matrix, var_threshold = 0.01) {
  cat("Filtering low variance genes...\n")

  gene_vars <- apply(gene_matrix, 2, var, na.rm = TRUE)
  high_var_genes <- gene_vars >= var_threshold

  cat(sprintf("Removed %d low variance genes (threshold: %.4f)\n",
              sum(!high_var_genes), var_threshold))
  cat(sprintf("Retained %d genes\n", sum(high_var_genes)))

  return(gene_matrix[, high_var_genes, drop = FALSE])
}


#' Handle Missing Values
#'
#' Imputes or removes missing values from gene expression data
#'
#' @param gene_matrix Matrix of gene expression
#' @param method Method: "remove", "mean", "median", "knn"
#' @param max_missing Maximum proportion of missing values per gene (for removal)
#' @return Processed gene matrix
handle_missing_values <- function(gene_matrix, method = "mean", max_missing = 0.2) {
  cat(sprintf("Handling missing values (method: %s)...\n", method))

  # Calculate missing proportion per gene
  missing_prop <- colMeans(is.na(gene_matrix))
  genes_to_keep <- missing_prop <= max_missing

  if (sum(!genes_to_keep) > 0) {
    cat(sprintf("Removing %d genes with >%.1f%% missing values\n",
                sum(!genes_to_keep), max_missing * 100))
    gene_matrix <- gene_matrix[, genes_to_keep, drop = FALSE]
  }

  # Count remaining missing values
  n_missing <- sum(is.na(gene_matrix))

  if (n_missing > 0) {
    cat(sprintf("Imputing %d missing values...\n", n_missing))

    if (method == "mean") {
      # Impute with column means
      for (j in 1:ncol(gene_matrix)) {
        gene_matrix[is.na(gene_matrix[, j]), j] <- mean(gene_matrix[, j], na.rm = TRUE)
      }
    } else if (method == "median") {
      # Impute with column medians
      for (j in 1:ncol(gene_matrix)) {
        gene_matrix[is.na(gene_matrix[, j]), j] <- median(gene_matrix[, j], na.rm = TRUE)
      }
    } else if (method == "knn") {
      if (requireNamespace("impute", quietly = TRUE)) {
        gene_matrix <- impute::impute.knn(as.matrix(gene_matrix))$data
      } else {
        warning("Package 'impute' not available. Using mean imputation instead.")
        gene_matrix <- handle_missing_values(gene_matrix, method = "mean", max_missing = max_missing)
      }
    } else if (method == "remove") {
      # Remove rows with any missing values
      complete_rows <- complete.cases(gene_matrix)
      cat(sprintf("Removing %d samples with missing values\n", sum(!complete_rows)))
      gene_matrix <- gene_matrix[complete_rows, , drop = FALSE]
    }
  } else {
    cat("No missing values found\n")
  }

  return(gene_matrix)
}


#' Normalize Gene Expression Data
#'
#' Applies normalization to gene expression data
#'
#' @param gene_matrix Matrix of gene expression (samples x genes)
#' @param method Normalization method: "zscore", "minmax", "quantile", "log2"
#' @return Normalized gene matrix
normalize_genes <- function(gene_matrix, method = "zscore") {
  cat(sprintf("Normalizing gene expression (method: %s)...\n", method))

  if (method == "zscore") {
    # Z-score normalization (standardization)
    gene_matrix <- scale(gene_matrix)

  } else if (method == "minmax") {
    # Min-max normalization to [0, 1]
    gene_matrix <- apply(gene_matrix, 2, function(x) {
      (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
    })

  } else if (method == "quantile") {
    # Quantile normalization
    if (requireNamespace("preprocessCore", quietly = TRUE)) {
      gene_matrix <- preprocessCore::normalize.quantiles(as.matrix(gene_matrix))
    } else {
      warning("Package 'preprocessCore' not available. Using z-score instead.")
      gene_matrix <- scale(gene_matrix)
    }

  } else if (method == "log2") {
    # Log2 transformation (add 1 to avoid log(0))
    gene_matrix <- log2(gene_matrix + 1)

  } else {
    warning(sprintf("Unknown normalization method: %s. Returning original data.", method))
  }

  return(gene_matrix)
}


#' Remove Outlier Samples
#'
#' Identifies and removes outlier samples based on distance metrics
#'
#' @param gene_matrix Matrix of gene expression (samples x genes)
#' @param method Method: "pca", "distance"
#' @param threshold Z-score threshold for outlier detection
#' @return List with cleaned matrix and outlier indices
remove_outliers <- function(gene_matrix, method = "pca", threshold = 3) {
  cat(sprintf("Detecting outliers (method: %s, threshold: %.1f)...\n", method, threshold))

  outliers <- c()

  if (method == "pca") {
    # Use PCA to detect outliers
    pca_result <- prcomp(gene_matrix, center = TRUE, scale. = TRUE)
    pc1_scores <- abs(scale(pca_result$x[, 1]))

    outliers <- which(pc1_scores > threshold)

  } else if (method == "distance") {
    # Use Euclidean distance from centroid
    centroid <- colMeans(gene_matrix)
    distances <- apply(gene_matrix, 1, function(x) sqrt(sum((x - centroid)^2)))
    z_distances <- abs(scale(distances))

    outliers <- which(z_distances > threshold)
  }

  if (length(outliers) > 0) {
    cat(sprintf("Identified %d outlier samples\n", length(outliers)))
    gene_matrix_clean <- gene_matrix[-outliers, , drop = FALSE]
  } else {
    cat("No outliers detected\n")
    gene_matrix_clean <- gene_matrix
  }

  result <- list(
    data = gene_matrix_clean,
    outlier_indices = outliers,
    n_outliers = length(outliers)
  )

  return(result)
}


#' Feature Selection by Correlation
#'
#' Removes highly correlated features to reduce multicollinearity
#'
#' @param gene_matrix Matrix of gene expression (samples x genes)
#' @param cor_threshold Correlation threshold (default: 0.9)
#' @return Filtered gene matrix
remove_correlated_features <- function(gene_matrix, cor_threshold = 0.9) {
  cat(sprintf("Removing highly correlated features (threshold: %.2f)...\n", cor_threshold))

  # Calculate correlation matrix
  cor_matrix <- cor(gene_matrix, use = "pairwise.complete.obs")

  # Find highly correlated pairs
  highly_corr <- which(abs(cor_matrix) > cor_threshold & upper.tri(cor_matrix), arr.ind = TRUE)

  if (nrow(highly_corr) > 0) {
    # Remove one gene from each highly correlated pair
    genes_to_remove <- unique(highly_corr[, 2])
    cat(sprintf("Removing %d highly correlated genes\n", length(genes_to_remove)))

    gene_matrix <- gene_matrix[, -genes_to_remove, drop = FALSE]
  } else {
    cat("No highly correlated features found\n")
  }

  return(gene_matrix)
}


#' Complete Preprocessing Pipeline
#'
#' Runs the full preprocessing pipeline on gene expression data
#'
#' @param gene_matrix Matrix of gene expression (samples x genes)
#' @param normalize Whether to normalize (default: TRUE)
#' @param norm_method Normalization method
#' @param filter_variance Whether to filter low variance genes
#' @param var_threshold Variance threshold
#' @param handle_na Whether to handle missing values
#' @param na_method Missing value method
#' @param remove_outliers_flag Whether to remove outliers
#' @param remove_corr Whether to remove correlated features
#' @param cor_threshold Correlation threshold
#' @return Preprocessed gene matrix
preprocess_pipeline <- function(gene_matrix,
                               normalize = TRUE,
                               norm_method = "zscore",
                               filter_variance = TRUE,
                               var_threshold = 0.01,
                               handle_na = TRUE,
                               na_method = "mean",
                               remove_outliers_flag = FALSE,
                               remove_corr = FALSE,
                               cor_threshold = 0.9) {

  cat("=== Starting Preprocessing Pipeline ===\n\n")

  original_dim <- dim(gene_matrix)
  cat(sprintf("Original dimensions: %d samples x %d genes\n\n", original_dim[1], original_dim[2]))

  # Step 1: Handle missing values
  if (handle_na) {
    gene_matrix <- handle_missing_values(gene_matrix, method = na_method)
    cat("\n")
  }

  # Step 2: Remove outlier samples
  if (remove_outliers_flag) {
    outlier_result <- remove_outliers(gene_matrix)
    gene_matrix <- outlier_result$data
    cat("\n")
  }

  # Step 3: Filter low variance genes
  if (filter_variance) {
    gene_matrix <- filter_low_variance(gene_matrix, var_threshold = var_threshold)
    cat("\n")
  }

  # Step 4: Remove highly correlated features
  if (remove_corr) {
    gene_matrix <- remove_correlated_features(gene_matrix, cor_threshold = cor_threshold)
    cat("\n")
  }

  # Step 5: Normalize
  if (normalize) {
    gene_matrix <- normalize_genes(gene_matrix, method = norm_method)
    cat("\n")
  }

  final_dim <- dim(gene_matrix)
  cat("=== Preprocessing Complete ===\n")
  cat(sprintf("Final dimensions: %d samples x %d genes\n", final_dim[1], final_dim[2]))
  cat(sprintf("Removed: %d samples, %d genes\n",
              original_dim[1] - final_dim[1],
              original_dim[2] - final_dim[2]))

  return(gene_matrix)
}


#' Split Data into Training and Testing Sets
#'
#' Splits data into training and testing sets
#'
#' @param x Feature matrix
#' @param y Response variable
#' @param train_ratio Proportion for training (default: 0.8)
#' @param seed Random seed for reproducibility
#' @return List with train and test data
train_test_split <- function(x, y, train_ratio = 0.8, seed = 123) {
  set.seed(seed)

  n <- nrow(x)
  train_size <- floor(train_ratio * n)
  train_indices <- sample(1:n, train_size)

  result <- list(
    x_train = x[train_indices, , drop = FALSE],
    x_test = x[-train_indices, , drop = FALSE],
    y_train = y[train_indices],
    y_test = y[-train_indices],
    train_indices = train_indices,
    test_indices = (1:n)[-train_indices]
  )

  cat(sprintf("Split data: %d training samples, %d testing samples\n",
              length(train_indices), n - length(train_indices)))

  return(result)
}
