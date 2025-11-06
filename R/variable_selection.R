library(glmnet)
library(randomForest)
library(caret)
library(MASS)

lasso_selection <- function(x, y, alpha = 1, nfolds = 10) {
  cat("Running LASSO variable selection...\n")

  # Perform cross-validation to find optimal lambda
  cv_fit <- cv.glmnet(x, y, alpha = alpha, nfolds = nfolds, family = "gaussian")

  # Extract coefficients at optimal lambda
  coef_lambda_min <- coef(cv_fit, s = "lambda.min")
  coef_lambda_1se <- coef(cv_fit, s = "lambda.1se")

  # Get selected variables (non-zero coefficients)
  selected_vars_min <- rownames(coef_lambda_min)[which(coef_lambda_min != 0)]
  selected_vars_1se <- rownames(coef_lambda_1se)[which(coef_lambda_1se != 0)]

  # Remove intercept
  selected_vars_min <- selected_vars_min[selected_vars_min != "(Intercept)"]
  selected_vars_1se <- selected_vars_1se[selected_vars_1se != "(Intercept)"]

  cat(sprintf("Variables selected (lambda.min): %d\n", length(selected_vars_min)))
  cat(sprintf("Variables selected (lambda.1se): %d\n", length(selected_vars_1se)))

  result <- list(
    cv_fit = cv_fit,
    selected_vars_min = selected_vars_min,
    selected_vars_1se = selected_vars_1se,
    lambda_min = cv_fit$lambda.min,
    lambda_1se = cv_fit$lambda.1se,
    coefficients_min = coef_lambda_min,
    coefficients_1se = coef_lambda_1se
  )

  return(result)
}

elastic_net_selection <- function(x, y, alpha = 0.5, nfolds = 10) {
  cat("Running Elastic Net variable selection...\n")

  # Use LASSO selection with specified alpha
  result <- lasso_selection(x, y, alpha = alpha, nfolds = nfolds)

  return(result)
}

rf_variable_importance <- function(x, y, ntree = 500, top_n = NULL, importance_threshold = NULL) {
  cat("Running Random Forest variable importance...\n")

  # Convert to data frame if matrix
  if (is.matrix(x)) {
    x <- as.data.frame(x)
  }

  # Fit random forest
  rf_model <- randomForest(x, y, ntree = ntree, importance = TRUE)

  # Get variable importance
  importance_scores <- importance(rf_model)
  importance_df <- data.frame(
    variable = rownames(importance_scores),
    importance = importance_scores[, "%IncMSE"],
    stringsAsFactors = FALSE
  )

  # Sort by importance
  importance_df <- importance_df[order(-importance_df$importance), ]

  # Select variables
  if (!is.null(top_n)) {
    selected_vars <- importance_df$variable[1:min(top_n, nrow(importance_df))]
    cat(sprintf("Selected top %d variables\n", length(selected_vars)))
  } else if (!is.null(importance_threshold)) {
    selected_vars <- importance_df$variable[importance_df$importance >= importance_threshold]
    cat(sprintf("Selected %d variables above threshold %.2f\n",
                length(selected_vars), importance_threshold))
  } else {
    selected_vars <- importance_df$variable
    cat(sprintf("Returning all %d variables with importance scores\n", length(selected_vars)))
  }

  result <- list(
    model = rf_model,
    importance = importance_df,
    selected_vars = selected_vars,
    top_importance = head(importance_df, 20)
  )

  return(result)
}

stepwise_selection <- function(x, y, direction = "both", max_vars = NULL) {
  cat(sprintf("Running stepwise variable selection (direction: %s)...\n", direction))

  # Create data frame
  data <- as.data.frame(cbind(y = y, x))

  # Start with null or full model depending on direction
  if (direction == "forward") {
    null_model <- lm(y ~ 1, data = data)
    full_formula <- as.formula(paste("y ~", paste(colnames(x), collapse = " + ")))

    step_model <- step(null_model,
                       scope = list(lower = null_model, upper = full_formula),
                       direction = "forward",
                       trace = 0)
  } else {
    full_formula <- as.formula(paste("y ~", paste(colnames(x), collapse = " + ")))
    full_model <- lm(full_formula, data = data)

    step_model <- step(full_model, direction = direction, trace = 0)
  }

  # Extract selected variables
  selected_vars <- names(coef(step_model))
  selected_vars <- selected_vars[selected_vars != "(Intercept)"]

  # Apply max_vars constraint if specified
  if (!is.null(max_vars) && length(selected_vars) > max_vars) {
    # Keep variables with largest absolute coefficients
    coefs <- coef(step_model)[-1]  # Remove intercept
    top_indices <- order(abs(coefs), decreasing = TRUE)[1:max_vars]
    selected_vars <- selected_vars[top_indices]
  }

  cat(sprintf("Selected %d variables\n", length(selected_vars)))

  result <- list(
    model = step_model,
    selected_vars = selected_vars,
    aic = AIC(step_model),
    bic = BIC(step_model)
  )

  return(result)
}

boruta_selection <- function(x, y, maxRuns = 100) {
  if (!requireNamespace("Boruta", quietly = TRUE)) {
    stop("Package 'Boruta' is required. Install with: install.packages('Boruta')")
  }

  cat("Running Boruta variable selection...\n")

  # Convert to data frame
  data <- as.data.frame(cbind(y = y, x))

  # Run Boruta
  boruta_result <- Boruta::Boruta(y ~ ., data = data, maxRuns = maxRuns)

  # Get confirmed and tentative variables
  confirmed_vars <- names(boruta_result$finalDecision[boruta_result$finalDecision == "Confirmed"])
  tentative_vars <- names(boruta_result$finalDecision[boruta_result$finalDecision == "Tentative"])

  cat(sprintf("Confirmed variables: %d\n", length(confirmed_vars)))
  cat(sprintf("Tentative variables: %d\n", length(tentative_vars)))

  result <- list(
    boruta_object = boruta_result,
    confirmed_vars = confirmed_vars,
    tentative_vars = tentative_vars,
    all_selected = c(confirmed_vars, tentative_vars)
  )

  return(result)
}

compare_methods <- function(x, y, methods = c("lasso", "elasticnet", "rf", "stepwise")) {
  cat("Comparing variable selection methods...\n\n")

  results <- list()

  if ("lasso" %in% methods) {
    results$lasso <- lasso_selection(x, y)
    cat("\n")
  }

  if ("elasticnet" %in% methods) {
    results$elasticnet <- elastic_net_selection(x, y, alpha = 0.5)
    cat("\n")
  }

  if ("rf" %in% methods) {
    results$rf <- rf_variable_importance(x, y, top_n = 50)
    cat("\n")
  }

  if ("stepwise" %in% methods) {
    # Use subset of variables for stepwise to avoid computational issues
    n_vars <- min(50, ncol(x))
    x_subset <- x[, 1:n_vars]
    results$stepwise <- stepwise_selection(x_subset, y)
    cat("\n")
  }

  # Summary
  cat("=== Method Comparison Summary ===\n")
  if ("lasso" %in% names(results)) {
    cat(sprintf("LASSO: %d variables\n", length(results$lasso$selected_vars_1se)))
  }
  if ("elasticnet" %in% names(results)) {
    cat(sprintf("Elastic Net: %d variables\n", length(results$elasticnet$selected_vars_1se)))
  }
  if ("rf" %in% names(results)) {
    cat(sprintf("Random Forest: %d variables\n", length(results$rf$selected_vars)))
  }
  if ("stepwise" %in% names(results)) {
    cat(sprintf("Stepwise: %d variables\n", length(results$stepwise$selected_vars)))
  }

  return(results)
}
