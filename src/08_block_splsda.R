# Block sPLS-DA Analysis Script
# Implements multiblock sparse Partial Least Squares Discriminant Analysis

suppressMessages({
  library(mixOmics)
  library(dplyr)
  library(ggplot2)
  library(RColorBrewer)
})

cat("Loading Block sPLS-DA analysis functions...\n")

# Load preprocessed data
multiomics_dir <- file.path(data_dir, "multiomics")
load(file.path(multiomics_dir, "train_processed.RData"))
load(file.path(multiomics_dir, "test_processed.RData"))

# Function to determine optimal number of components
tune_block_splsda_components <- function(X, Y, design, validation = "Mfold", 
                                       folds = 5, nrepeat = 10, max_comp = 5) {
  cat("Tuning number of components for Block sPLS-DA...\n")
  
  # Perform cross-validation to select optimal number of components
  tune_result <- tune.block.splsda(
    X = X,
    Y = Y,
    ncomp = max_comp,
    design = design,
    validation = validation,
    folds = folds,
    nrepeat = nrepeat,
    measure = "BER"  # Balanced Error Rate
  )
  
  # Extract optimal number of components
  optimal_ncomp <- tune_result$choice.ncomp$ncomp
  
  cat("Optimal number of components:", optimal_ncomp, "\n")
  
  # Save tuning results
  save(tune_result, file = file.path(multiomics_dir, "component_tuning.RData"))
  
  return(list(
    tune_result = tune_result,
    optimal_ncomp = optimal_ncomp
  ))
}

# Function to tune keepX parameters (variable selection)
tune_keepx_parameters <- function(X, Y, design, ncomp, validation = "Mfold", 
                                folds = 5, nrepeat = 10) {
  cat("Tuning keepX parameters for variable selection...\n")
  
  # Define test values for keepX
  test_keepx <- list()
  for (block_name in names(X)) {
    if (block_name == "mRNA") {
      test_keepx[[block_name]] <- c(5, 10, 15, 20, 25, 30)
    } else if (block_name == "miRNA") {
      test_keepx[[block_name]] <- c(5, 10, 15, 20, 25)
    } else if (block_name == "protein") {
      test_keepx[[block_name]] <- c(5, 10, 15, 20)
    }
  }
  
  # Perform tuning
  tune_keepx_result <- tune.block.splsda(
    X = X,
    Y = Y,
    ncomp = ncomp,
    design = design,
    validation = validation,
    folds = folds,
    nrepeat = nrepeat,
    test.keepX = test_keepx,
    measure = "BER"
  )
  
  # Extract optimal keepX values
  optimal_keepx <- tune_keepx_result$choice.keepX
  
  cat("Optimal keepX parameters:\n")
  for (block_name in names(optimal_keepx)) {
    cat("  ", block_name, ": ", paste(optimal_keepx[[block_name]], collapse = ", "), "\n")
  }
  
  # Save tuning results
  save(tune_keepx_result, file = file.path(multiomics_dir, "keepx_tuning.RData"))
  
  return(list(
    tune_result = tune_keepx_result,
    optimal_keepx = optimal_keepx
  ))
}

# Function to build final Block sPLS-DA model
build_block_splsda_model <- function(X, Y, design, ncomp, keepX) {
  cat("Building final Block sPLS-DA model...\n")
  
  # Build the model
  block_splsda_model <- block.splsda(
    X = X,
    Y = Y,
    ncomp = ncomp,
    design = design,
    keepX = keepX
  )
  
  # Extract model information
  model_info <- list(
    ncomp = ncomp,
    keepX = keepX,
    design = design,
    explained_variance = block_splsda_model$explained_variance,
    loadings = block_splsda_model$loadings,
    variates = block_splsda_model$variates
  )
  
  cat("Model built with", ncomp, "components\n")
  cat("Explained variance per component:\n")
  for (i in 1:ncomp) {
    cat("  Component", i, ":\n")
    for (block_name in names(X)) {
      if (block_name %in% names(block_splsda_model$explained_variance)) {
        var_explained <- block_splsda_model$explained_variance[[block_name]][i]
        cat("    ", block_name, ": ", round(var_explained * 100, 2), "%\n")
      }
    }
  }
  
  return(list(
    model = block_splsda_model,
    info = model_info
  ))
}

# Function to evaluate model performance
evaluate_model_performance <- function(model, X_test, Y_test, design) {
  cat("Evaluating model performance...\n")
  
  # Make predictions on test data
  if (!is.null(X_test) && !is.null(Y_test)) {
    predictions <- predict(model, newdata = X_test, design = design)
    
    # Calculate prediction accuracy
    predicted_classes <- predictions$WeightedVote$centroids.dist[, 1]
    actual_classes <- Y_test
    
    # Calculate accuracy
    accuracy <- sum(predicted_classes == actual_classes) / length(actual_classes)
    
    # Calculate confusion matrix
    confusion_matrix <- table(Predicted = predicted_classes, Actual = actual_classes)
    
    # Calculate per-class metrics
    class_metrics <- list()
    for (class_name in levels(actual_classes)) {
      tp <- confusion_matrix[class_name, class_name]
      fp <- sum(confusion_matrix[class_name, ]) - tp
      fn <- sum(confusion_matrix[, class_name]) - tp
      tn <- sum(confusion_matrix) - tp - fp - fn
      
      precision <- tp / (tp + fp)
      recall <- tp / (tp + fn)
      f1 <- 2 * (precision * recall) / (precision + recall)
      
      class_metrics[[class_name]] <- list(
        precision = precision,
        recall = recall,
        f1 = f1
      )
    }
    
    performance_results <- list(
      accuracy = accuracy,
      confusion_matrix = confusion_matrix,
      class_metrics = class_metrics,
      predictions = predictions
    )
    
    cat("Model accuracy:", round(accuracy * 100, 2), "%\n")
    cat("Confusion matrix:\n")
    print(confusion_matrix)
    
  } else {
    cat("Test data not available for performance evaluation\n")
    performance_results <- NULL
  }
  
  return(performance_results)
}

# Function to extract important variables
extract_important_variables <- function(model, top_n = 20) {
  cat("Extracting important variables...\n")
  
  important_vars <- list()
  
  for (block_name in names(model$model$X)) {
    block_loadings <- model$model$loadings[[block_name]]
    
    # Calculate variable importance scores
    var_importance <- apply(abs(block_loadings), 1, sum)
    
    # Get top variables
    top_vars <- names(sort(var_importance, decreasing = TRUE))[1:min(top_n, length(var_importance))]
    
    important_vars[[block_name]] <- list(
      variables = top_vars,
      importance_scores = var_importance[top_vars],
      loadings = block_loadings[top_vars, , drop = FALSE]
    )
    
    cat("Top", length(top_vars), "variables in", block_name, ":\n")
    for (i in 1:min(10, length(top_vars))) {
      cat("  ", top_vars[i], ": ", round(var_importance[top_vars[i]], 3), "\n")
    }
  }
  
  return(important_vars)
}

# Function to perform cross-validation assessment
perform_cv_assessment <- function(X, Y, design, ncomp, keepX, folds = 5, nrepeat = 10) {
  cat("Performing cross-validation assessment...\n")
  
  # Perform repeated cross-validation
  cv_results <- perf(
    block.splsda(X = X, Y = Y, ncomp = ncomp, design = design, keepX = keepX),
    validation = "Mfold",
    folds = folds,
    nrepeat = nrepeat,
    measure = c("BER", "overall")
  )
  
  # Extract performance metrics
  cv_performance <- list(
    BER = cv_results$error.rate$BER,
    overall = cv_results$error.rate$overall,
    error_rate_per_class = cv_results$error.rate.class
  )
  
  cat("Cross-validation results:\n")
  cat("  BER (Balanced Error Rate):", round(cv_performance$BER[ncomp], 3), "\n")
  cat("  Overall error rate:", round(cv_performance$overall[ncomp], 3), "\n")
  
  return(cv_performance)
}

# Main analysis workflow
cat("=== Block sPLS-DA Analysis ===\n")

# Extract data for analysis
X_train <- train_processed$X
Y_train <- train_processed$Y
design_matrix <- train_processed$design

X_test <- test_processed$X
Y_test <- test_processed$Y

# Step 1: Tune number of components
component_tuning <- tune_block_splsda_components(
  X = X_train,
  Y = Y_train,
  design = design_matrix,
  max_comp = 3
)

# Step 2: Tune keepX parameters
keepx_tuning <- tune_keepx_parameters(
  X = X_train,
  Y = Y_train,
  design = design_matrix,
  ncomp = component_tuning$optimal_ncomp
)

# Step 3: Build final model
final_model <- build_block_splsda_model(
  X = X_train,
  Y = Y_train,
  design = design_matrix,
  ncomp = component_tuning$optimal_ncomp,
  keepX = keepx_tuning$optimal_keepx
)

# Step 4: Evaluate performance
performance_results <- evaluate_model_performance(
  model = final_model$model,
  X_test = X_test,
  Y_test = Y_test,
  design = design_matrix
)

# Step 5: Extract important variables
important_variables <- extract_important_variables(final_model, top_n = 30)

# Step 6: Cross-validation assessment
cv_assessment <- perform_cv_assessment(
  X = X_train,
  Y = Y_train,
  design = design_matrix,
  ncomp = component_tuning$optimal_ncomp,
  keepX = keepx_tuning$optimal_keepx
)

# Compile all results
block_splsda_results <- list(
  component_tuning = component_tuning,
  keepx_tuning = keepx_tuning,
  final_model = final_model,
  performance = performance_results,
  important_variables = important_variables,
  cv_assessment = cv_assessment
)

# Save all results
save(block_splsda_results, file = file.path(multiomics_dir, "block_splsda_results.RData"))

# Create summary tables
create_results_summary <- function(results) {
  # Model summary
  model_summary <- data.frame(
    Component = 1:results$final_model$info$ncomp,
    stringsAsFactors = FALSE
  )
  
  # Add explained variance for each block
  for (block_name in names(results$final_model$model$X)) {
    if (block_name %in% names(results$final_model$model$explained_variance)) {
      var_col <- paste0(block_name, "_variance")
      model_summary[[var_col]] <- results$final_model$model$explained_variance[[block_name]]
    }
  }
  
  # Performance summary
  if (!is.null(results$performance)) {
    performance_summary <- data.frame(
      Metric = c("Accuracy", "BER"),
      Value = c(results$performance$accuracy, results$cv_assessment$BER[results$final_model$info$ncomp])
    )
  } else {
    performance_summary <- NULL
  }
  
  # Important variables summary
  important_vars_summary <- list()
  for (block_name in names(results$important_variables)) {
    vars_df <- data.frame(
      Variable = results$important_variables[[block_name]]$variables,
      Importance = results$important_variables[[block_name]]$importance_scores,
      Block = block_name
    )
    important_vars_summary[[block_name]] <- vars_df
  }
  
  return(list(
    model_summary = model_summary,
    performance_summary = performance_summary,
    important_vars_summary = important_vars_summary
  ))
}

# Generate summaries
results_summary <- create_results_summary(block_splsda_results)

# Save summary tables
write.csv(results_summary$model_summary, 
          file = file.path(tables_dir, "block_splsda_model_summary.csv"), 
          row.names = FALSE)

if (!is.null(results_summary$performance_summary)) {
  write.csv(results_summary$performance_summary, 
            file = file.path(tables_dir, "block_splsda_performance.csv"), 
            row.names = FALSE)
}

for (block_name in names(results_summary$important_vars_summary)) {
  write.csv(results_summary$important_vars_summary[[block_name]], 
            file = file.path(tables_dir, paste0("block_splsda_important_vars_", block_name, ".csv")), 
            row.names = FALSE)
}

cat("Block sPLS-DA analysis completed successfully!\n")
cat("Results saved to:", multiomics_dir, "\n")
cat("Summary tables saved to:", tables_dir, "\n")