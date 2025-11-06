library(ggplot2)
library(pheatmap)

plot_variable_importance <- function(importance_df, top_n = 20, title = "Variable Importance") {
  # Select top N variables
  top_vars <- head(importance_df, top_n)

  # Create plot
  p <- ggplot(top_vars, aes(x = reorder(variable, importance), y = importance)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    labs(title = title,
         x = "Variables (Genes)",
         y = "Importance Score") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 8))

  return(p)
}

plot_lasso_path <- function(lasso_result) {
  cv_fit <- lasso_result$cv_fit

  par(mfrow = c(1, 2))

  # Plot coefficient paths
  plot(cv_fit$glmnet.fit, xvar = "lambda", label = TRUE, main = "LASSO Coefficient Paths")
  abline(v = log(cv_fit$lambda.min), lty = 2, col = "red")
  abline(v = log(cv_fit$lambda.1se), lty = 2, col = "blue")

  # Plot cross-validation curve
  plot(cv_fit, main = "Cross-Validation Curve")

  par(mfrow = c(1, 1))
}

plot_gene_heatmap <- function(gene_matrix, selected_genes, annotation_col = NULL,
                              title = "Selected Gene Expression") {

  # Subset to selected genes
  if (is.matrix(gene_matrix)) {
    selected_data <- gene_matrix[, selected_genes, drop = FALSE]
  } else {
    selected_data <- as.matrix(gene_matrix[, selected_genes, drop = FALSE])
  }

  # Create heatmap
  pheatmap(t(selected_data),
           scale = "row",
           clustering_distance_rows = "euclidean",
           clustering_distance_cols = "euclidean",
           annotation_col = annotation_col,
           main = title,
           fontsize_row = 8,
           fontsize_col = 8)
}

calculate_performance <- function(y_true, y_pred) {
  # Remove NA values
  valid_idx <- !is.na(y_true) & !is.na(y_pred)
  y_true <- y_true[valid_idx]
  y_pred <- y_pred[valid_idx]

  # Calculate metrics
  mse <- mean((y_true - y_pred)^2)
  rmse <- sqrt(mse)
  mae <- mean(abs(y_true - y_pred))

  # R-squared
  ss_res <- sum((y_true - y_pred)^2)
  ss_tot <- sum((y_true - mean(y_true))^2)
  r_squared <- 1 - (ss_res / ss_tot)

  # Correlation
  correlation <- cor(y_true, y_pred)

  metrics <- list(
    MSE = mse,
    RMSE = rmse,
    MAE = mae,
    R_squared = r_squared,
    Correlation = correlation,
    n = length(y_true)
  )

  return(metrics)
}

print_performance <- function(metrics, method_name = NULL) {
  if (!is.null(method_name)) {
    cat(sprintf("=== Performance Metrics: %s ===\n", method_name))
  } else {
    cat("=== Performance Metrics ===\n")
  }

  cat(sprintf("Sample size: %d\n", metrics$n))
  cat(sprintf("RMSE:        %.4f\n", metrics$RMSE))
  cat(sprintf("MAE:         %.4f\n", metrics$MAE))
  cat(sprintf("R-squared:   %.4f\n", metrics$R_squared))
  cat(sprintf("Correlation: %.4f\n", metrics$Correlation))
  cat("\n")
}

validate_selection <- function(x_train, y_train, x_test, y_test,
                               selected_vars, method = "lm") {

  cat(sprintf("Validating %d selected variables using %s...\n", length(selected_vars), method))

  # Subset to selected variables
  x_train_subset <- x_train[, selected_vars, drop = FALSE]
  x_test_subset <- x_test[, selected_vars, drop = FALSE]

  # Train model
  if (method == "lm") {
    data_train <- as.data.frame(cbind(y = y_train, x_train_subset))
    model <- lm(y ~ ., data = data_train)

    # Predict
    data_test <- as.data.frame(x_test_subset)
    y_pred <- predict(model, newdata = data_test)

  } else if (method == "glmnet") {
    model <- glmnet::glmnet(x_train_subset, y_train, alpha = 1)
    y_pred <- predict(model, newx = x_test_subset, s = 0.01)

  } else if (method == "rf") {
    data_train <- as.data.frame(cbind(y = y_train, x_train_subset))
    model <- randomForest::randomForest(y ~ ., data = data_train)

    data_test <- as.data.frame(x_test_subset)
    y_pred <- predict(model, newdata = data_test)
  }

  # Calculate performance
  metrics <- calculate_performance(y_test, y_pred)

  result <- list(
    model = model,
    predictions = y_pred,
    metrics = metrics
  )

  return(result)
}

compare_performance <- function(comparison_results, x_train, y_train, x_test, y_test) {
  cat("=== Comparing Method Performance ===\n\n")

  performance_df <- data.frame(
    Method = character(),
    N_Variables = integer(),
    RMSE = numeric(),
    MAE = numeric(),
    R_squared = numeric(),
    stringsAsFactors = FALSE
  )

  for (method_name in names(comparison_results)) {
    result <- comparison_results[[method_name]]

    # Get selected variables
    if (method_name %in% c("lasso", "elasticnet")) {
      selected_vars <- result$selected_vars_1se
    } else if (method_name == "rf") {
      selected_vars <- result$selected_vars
    } else if (method_name == "stepwise") {
      selected_vars <- result$selected_vars
    }

    # Skip if no variables selected
    if (length(selected_vars) == 0) {
      cat(sprintf("Skipping %s: no variables selected\n", method_name))
      next
    }

    # Validate
    validation <- validate_selection(x_train, y_train, x_test, y_test,
                                     selected_vars, method = "lm")

    print_performance(validation$metrics, method_name)

    # Add to results
    performance_df <- rbind(performance_df, data.frame(
      Method = method_name,
      N_Variables = length(selected_vars),
      RMSE = validation$metrics$RMSE,
      MAE = validation$metrics$MAE,
      R_squared = validation$metrics$R_squared,
      stringsAsFactors = FALSE
    ))
  }

  return(performance_df)
}

export_selected_vars <- function(selected_vars, output_file, method_name = NULL) {
  cat(sprintf("Exporting %d selected variables to %s\n", length(selected_vars), output_file))

  # Create data frame
  export_df <- data.frame(
    Variable = selected_vars,
    Index = 1:length(selected_vars),
    stringsAsFactors = FALSE
  )

  if (!is.null(method_name)) {
    export_df$Method <- method_name
  }

  # Write to file
  write.csv(export_df, output_file, row.names = FALSE)

  cat("Export complete\n")
}

generate_report <- function(results, method_name, output_file = NULL) {
  report <- c()
  report <- c(report, paste0("=== Variable Selection Report: ", method_name, " ===\n"))
  report <- c(report, paste0("Date: ", Sys.time(), "\n\n"))

  if (method_name %in% c("lasso", "elasticnet")) {
    report <- c(report, sprintf("Lambda (min): %.6f\n", results$lambda_min))
    report <- c(report, sprintf("Lambda (1se): %.6f\n", results$lambda_1se))
    report <- c(report, sprintf("Variables selected (lambda.min): %d\n", length(results$selected_vars_min)))
    report <- c(report, sprintf("Variables selected (lambda.1se): %d\n\n", length(results$selected_vars_1se)))

    report <- c(report, "Selected variables (lambda.1se):\n")
    for (var in results$selected_vars_1se) {
      report <- c(report, sprintf("  - %s\n", var))
    }

  } else if (method_name == "rf") {
    report <- c(report, sprintf("Number of trees: %d\n", results$model$ntree))
    report <- c(report, sprintf("Variables selected: %d\n\n", length(results$selected_vars)))

    report <- c(report, "Top 20 variables by importance:\n")
    top_vars <- head(results$importance, 20)
    for (i in 1:nrow(top_vars)) {
      report <- c(report, sprintf("  %2d. %s (%.4f)\n",
                                 i, top_vars$variable[i], top_vars$importance[i]))
    }
  }

  # Print to console
  cat(paste(report, collapse = ""))

  # Optionally write to file
  if (!is.null(output_file)) {
    writeLines(report, output_file)
    cat(sprintf("\nReport saved to: %s\n", output_file))
  }

  invisible(report)
}

create_enrichment_input <- function(selected_genes, output_file) {
  cat(sprintf("Creating gene set enrichment input with %d genes\n", length(selected_genes)))

  writeLines(selected_genes, output_file)

  cat(sprintf("File created: %s\n", output_file))
  cat("This file can be used with tools like DAVID, Enrichr, or g:Profiler\n")
}
