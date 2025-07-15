# DIABLO Analysis Script
# Data Integration Analysis for Biomarker discovery using Latent variable approaches

suppressMessages({
  library(mixOmics)
  library(dplyr)
  library(ggplot2)
  library(RColorBrewer)
  library(gridExtra)
})

cat("Loading DIABLO analysis functions...\n")

# Load preprocessed data
multiomics_dir <- file.path(data_dir, "multiomics")
load(file.path(multiomics_dir, "train_processed.RData"))
load(file.path(multiomics_dir, "test_processed.RData"))

# Function to run DIABLO analysis
run_diablo_analysis <- function(X, Y, design, ncomp = 2, validation = "Mfold", 
                               folds = 5, nrepeat = 10) {
  cat("Running DIABLO analysis...\n")
  
  # Step 1: Tune number of components
  cat("  Step 1: Tuning number of components...\n")
  diablo_tune <- tune(
    X = X,
    Y = Y,
    ncomp = ncomp,
    design = design,
    validation = validation,
    folds = folds,
    nrepeat = nrepeat,
    measure = "BER"
  )
  
  optimal_ncomp <- diablo_tune$choice.ncomp$ncomp
  cat("    Optimal number of components:", optimal_ncomp, "\n")
  
  # Step 2: Tune keepX parameters
  cat("  Step 2: Tuning keepX parameters...\n")
  
  # Define test keepX values
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
  
  diablo_keepx <- tune.block.splsda(
    X = X,
    Y = Y,
    ncomp = optimal_ncomp,
    design = design,
    validation = validation,
    folds = folds,
    nrepeat = nrepeat,
    test.keepX = test_keepx,
    measure = "BER"
  )
  
  optimal_keepx <- diablo_keepx$choice.keepX
  cat("    Optimal keepX values:\n")
  for (block_name in names(optimal_keepx)) {
    cat("      ", block_name, ": ", paste(optimal_keepx[[block_name]], collapse = ", "), "\n")
  }
  
  # Step 3: Build final DIABLO model
  cat("  Step 3: Building final DIABLO model...\n")
  diablo_model <- block.splsda(
    X = X,
    Y = Y,
    ncomp = optimal_ncomp,
    design = design,
    keepX = optimal_keepx
  )
  
  # Step 4: Evaluate model performance
  cat("  Step 4: Evaluating model performance...\n")
  diablo_perf <- perf(
    diablo_model,
    validation = validation,
    folds = folds,
    nrepeat = nrepeat,
    measure = c("BER", "overall")
  )
  
  cat("    Model performance:\n")
  cat("      BER:", round(diablo_perf$error.rate$BER[optimal_ncomp], 3), "\n")
  cat("      Overall error rate:", round(diablo_perf$error.rate$overall[optimal_ncomp], 3), "\n")
  
  return(list(
    component_tuning = diablo_tune,
    keepx_tuning = diablo_keepx,
    model = diablo_model,
    performance = diablo_perf,
    optimal_ncomp = optimal_ncomp,
    optimal_keepx = optimal_keepx
  ))
}

# Function to make predictions with DIABLO
predict_diablo <- function(diablo_model, X_test, Y_test = NULL, design) {
  cat("Making predictions with DIABLO model...\n")
  
  # Handle missing data blocks
  X_test_filtered <- list()
  for (block_name in names(diablo_model$X)) {
    if (block_name %in% names(X_test) && !is.null(X_test[[block_name]])) {
      X_test_filtered[[block_name]] <- X_test[[block_name]]
    } else {
      cat("  Missing data block:", block_name, "\n")
    }
  }
  
  if (length(X_test_filtered) == 0) {
    cat("  No valid data blocks for prediction\n")
    return(NULL)
  }
  
  # Make predictions
  predictions <- predict(diablo_model, newdata = X_test_filtered, design = design)
  
  # Extract different prediction methods
  prediction_results <- list(
    centroids_dist = predictions$centroids.dist,
    majority_vote = predictions$MajorityVote,
    weighted_vote = predictions$WeightedVote
  )
  
  # Evaluate predictions if true labels are available
  if (!is.null(Y_test)) {
    cat("  Evaluating prediction performance...\n")
    
    # Get predictions for each method
    pred_centroids <- predictions$centroids.dist[, 1]
    pred_majority <- predictions$MajorityVote[, 1]
    pred_weighted <- predictions$WeightedVote[, 1]
    
    # Calculate accuracies
    acc_centroids <- sum(pred_centroids == Y_test) / length(Y_test)
    acc_majority <- sum(pred_majority == Y_test) / length(Y_test)
    acc_weighted <- sum(pred_weighted == Y_test) / length(Y_test)
    
    cat("    Prediction accuracies:\n")
    cat("      Centroids distance:", round(acc_centroids * 100, 2), "%\n")
    cat("      Majority vote:", round(acc_majority * 100, 2), "%\n")
    cat("      Weighted vote:", round(acc_weighted * 100, 2), "%\n")
    
    # Create confusion matrices
    conf_centroids <- table(Predicted = pred_centroids, Actual = Y_test)
    conf_majority <- table(Predicted = pred_majority, Actual = Y_test)
    conf_weighted <- table(Predicted = pred_weighted, Actual = Y_test)
    
    prediction_results$evaluation <- list(
      accuracy_centroids = acc_centroids,
      accuracy_majority = acc_majority,
      accuracy_weighted = acc_weighted,
      confusion_centroids = conf_centroids,
      confusion_majority = conf_majority,
      confusion_weighted = conf_weighted
    )
  }
  
  return(prediction_results)
}

# Function to extract biomarker signatures
extract_biomarker_signatures <- function(diablo_model, top_n = 20) {
  cat("Extracting biomarker signatures...\n")
  
  biomarkers <- list()
  
  for (block_name in names(diablo_model$X)) {
    block_loadings <- diablo_model$loadings[[block_name]]
    
    # Calculate stability/importance scores
    stability_scores <- apply(abs(block_loadings), 1, sum)
    
    # Get top biomarkers
    top_biomarkers <- names(sort(stability_scores, decreasing = TRUE))[1:min(top_n, length(stability_scores))]
    
    biomarkers[[block_name]] <- list(
      biomarkers = top_biomarkers,
      stability_scores = stability_scores[top_biomarkers],
      loadings = block_loadings[top_biomarkers, , drop = FALSE]
    )
    
    cat("  Top", length(top_biomarkers), "biomarkers in", block_name, ":\n")
    for (i in 1:min(10, length(top_biomarkers))) {
      cat("    ", top_biomarkers[i], ": ", round(stability_scores[top_biomarkers[i]], 3), "\n")
    }
  }
  
  return(biomarkers)
}

# Function to perform feature correlation analysis
analyze_feature_correlations <- function(diablo_model, min_cor = 0.7) {
  cat("Analyzing feature correlations across data blocks...\n")
  
  # Get variates (component scores)
  variates <- diablo_model$variates
  
  # Calculate correlations between variates of different blocks
  block_correlations <- list()
  block_names <- names(variates)
  
  for (i in 1:(length(block_names) - 1)) {
    for (j in (i + 1):length(block_names)) {
      block1 <- block_names[i]
      block2 <- block_names[j]
      
      # Calculate correlation matrix
      cor_matrix <- cor(variates[[block1]], variates[[block2]])
      
      block_correlations[[paste(block1, block2, sep = "_")]] <- cor_matrix
      
      cat("  Correlation between", block1, "and", block2, "variates:\n")
      print(round(cor_matrix, 3))
    }
  }
  
  # Find highly correlated features
  highly_correlated <- list()
  for (block_name in block_names) {
    block_loadings <- diablo_model$loadings[[block_name]]
    
    # Calculate pairwise correlations of loadings
    loading_cors <- cor(t(block_loadings))
    
    # Find highly correlated pairs
    high_cor_pairs <- which(abs(loading_cors) > min_cor & loading_cors != 1, arr.ind = TRUE)
    
    if (nrow(high_cor_pairs) > 0) {
      cor_pairs <- data.frame(
        feature1 = rownames(loading_cors)[high_cor_pairs[, 1]],
        feature2 = rownames(loading_cors)[high_cor_pairs[, 2]],
        correlation = loading_cors[high_cor_pairs],
        block = block_name
      )
      highly_correlated[[block_name]] <- cor_pairs
    }
  }
  
  return(list(
    block_correlations = block_correlations,
    highly_correlated_features = highly_correlated
  ))
}

# Function to assess model stability
assess_model_stability <- function(X, Y, design, ncomp, keepX, n_bootstrap = 100) {
  cat("Assessing model stability through bootstrap...\n")
  
  n_samples <- nrow(X[[1]])
  stability_results <- list()
  
  # Initialize stability matrices
  for (block_name in names(X)) {
    stability_results[[block_name]] <- matrix(0, 
                                             nrow = ncol(X[[block_name]]), 
                                             ncol = ncomp)
    rownames(stability_results[[block_name]]) <- colnames(X[[block_name]])
  }
  
  # Bootstrap sampling
  for (i in 1:n_bootstrap) {
    if (i %% 20 == 0) cat("  Bootstrap iteration", i, "/", n_bootstrap, "\n")
    
    # Sample with replacement
    boot_indices <- sample(1:n_samples, n_samples, replace = TRUE)
    
    # Create bootstrap datasets
    X_boot <- list()
    for (block_name in names(X)) {
      X_boot[[block_name]] <- X[[block_name]][boot_indices, ]
    }
    Y_boot <- Y[boot_indices]
    
    # Fit model on bootstrap sample
    tryCatch({
      boot_model <- block.splsda(
        X = X_boot,
        Y = Y_boot,
        ncomp = ncomp,
        design = design,
        keepX = keepX
      )
      
      # Update stability scores
      for (block_name in names(X)) {
        selected_features <- rownames(boot_model$loadings[[block_name]])
        selected_features <- selected_features[selected_features != ""]
        
        for (comp in 1:ncomp) {
          stability_results[[block_name]][selected_features, comp] <- 
            stability_results[[block_name]][selected_features, comp] + 1
        }
      }
    }, error = function(e) {
      cat("    Bootstrap iteration", i, "failed:", e$message, "\n")
    })
  }
  
  # Calculate stability percentages
  for (block_name in names(stability_results)) {
    stability_results[[block_name]] <- stability_results[[block_name]] / n_bootstrap
  }
  
  cat("  Stability assessment completed\n")
  return(stability_results)
}

# Main DIABLO analysis
cat("=== DIABLO Analysis ===\n")

# Extract data
X_train <- train_processed$X
Y_train <- train_processed$Y
design_matrix <- train_processed$design

X_test <- test_processed$X
Y_test <- test_processed$Y

# Run DIABLO analysis
diablo_results <- run_diablo_analysis(
  X = X_train,
  Y = Y_train,
  design = design_matrix,
  ncomp = 3
)

# Make predictions on test data
predictions <- predict_diablo(
  diablo_model = diablo_results$model,
  X_test = X_test,
  Y_test = Y_test,
  design = design_matrix
)

# Extract biomarker signatures
biomarkers <- extract_biomarker_signatures(diablo_results$model, top_n = 30)

# Analyze feature correlations
correlations <- analyze_feature_correlations(diablo_results$model)

# Assess model stability
stability <- assess_model_stability(
  X = X_train,
  Y = Y_train,
  design = design_matrix,
  ncomp = diablo_results$optimal_ncomp,
  keepX = diablo_results$optimal_keepx,
  n_bootstrap = 50
)

# Compile all results
diablo_final_results <- list(
  main_results = diablo_results,
  predictions = predictions,
  biomarkers = biomarkers,
  correlations = correlations,
  stability = stability
)

# Save results
save(diablo_final_results, file = file.path(multiomics_dir, "diablo_results.RData"))

# Create summary tables
create_diablo_summary_tables <- function(results) {
  # Performance summary
  perf_summary <- data.frame(
    Metric = c("BER", "Overall Error", "Centroids Accuracy", "Majority Vote Accuracy", "Weighted Vote Accuracy"),
    Value = c(
      results$main_results$performance$error.rate$BER[results$main_results$optimal_ncomp],
      results$main_results$performance$error.rate$overall[results$main_results$optimal_ncomp],
      if (!is.null(results$predictions$evaluation)) results$predictions$evaluation$accuracy_centroids else NA,
      if (!is.null(results$predictions$evaluation)) results$predictions$evaluation$accuracy_majority else NA,
      if (!is.null(results$predictions$evaluation)) results$predictions$evaluation$accuracy_weighted else NA
    )
  )
  
  # Biomarker summary
  biomarker_summary <- list()
  for (block_name in names(results$biomarkers)) {
    biomarker_df <- data.frame(
      Biomarker = results$biomarkers[[block_name]]$biomarkers,
      Stability_Score = results$biomarkers[[block_name]]$stability_scores,
      Block = block_name
    )
    biomarker_summary[[block_name]] <- biomarker_df
  }
  
  # Model configuration summary
  config_summary <- data.frame(
    Parameter = c("Components", paste0("keepX_", names(results$main_results$optimal_keepx))),
    Value = c(results$main_results$optimal_ncomp, unlist(results$main_results$optimal_keepx))
  )
  
  return(list(
    performance = perf_summary,
    biomarkers = biomarker_summary,
    configuration = config_summary
  ))
}

# Generate summary tables
summary_tables <- create_diablo_summary_tables(diablo_final_results)

# Save summary tables
write.csv(summary_tables$performance, 
          file = file.path(tables_dir, "diablo_performance_summary.csv"), 
          row.names = FALSE)

write.csv(summary_tables$configuration, 
          file = file.path(tables_dir, "diablo_configuration.csv"), 
          row.names = FALSE)

for (block_name in names(summary_tables$biomarkers)) {
  write.csv(summary_tables$biomarkers[[block_name]], 
            file = file.path(tables_dir, paste0("diablo_biomarkers_", block_name, ".csv")), 
            row.names = FALSE)
}

# Save prediction results if available
if (!is.null(predictions) && !is.null(predictions$evaluation)) {
  # Save confusion matrices
  write.csv(predictions$evaluation$confusion_centroids, 
            file = file.path(tables_dir, "diablo_confusion_centroids.csv"))
  write.csv(predictions$evaluation$confusion_majority, 
            file = file.path(tables_dir, "diablo_confusion_majority.csv"))
  write.csv(predictions$evaluation$confusion_weighted, 
            file = file.path(tables_dir, "diablo_confusion_weighted.csv"))
}

cat("DIABLO analysis completed successfully!\n")
cat("Results saved to:", multiomics_dir, "\n")
cat("Summary tables saved to:", tables_dir, "\n")

# Print final summary
cat("\n=== DIABLO Analysis Summary ===\n")
cat("Model configuration:\n")
cat("  Components:", diablo_final_results$main_results$optimal_ncomp, "\n")
cat("  keepX parameters:\n")
for (block_name in names(diablo_final_results$main_results$optimal_keepx)) {
  cat("    ", block_name, ":", paste(diablo_final_results$main_results$optimal_keepx[[block_name]], collapse = ", "), "\n")
}

cat("\nModel performance:\n")
cat("  Cross-validation BER:", round(diablo_final_results$main_results$performance$error.rate$BER[diablo_final_results$main_results$optimal_ncomp], 3), "\n")

if (!is.null(predictions) && !is.null(predictions$evaluation)) {
  cat("\nTest set predictions:\n")
  cat("  Centroids distance accuracy:", round(predictions$evaluation$accuracy_centroids * 100, 2), "%\n")
  cat("  Majority vote accuracy:", round(predictions$evaluation$accuracy_majority * 100, 2), "%\n")
  cat("  Weighted vote accuracy:", round(predictions$evaluation$accuracy_weighted * 100, 2), "%\n")
}

cat("\nTop biomarkers per block:\n")
for (block_name in names(biomarkers)) {
  cat("  ", block_name, ":", paste(biomarkers[[block_name]]$biomarkers[1:5], collapse = ", "), "\n")
}