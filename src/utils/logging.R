# Centralized Logging Utilities for TCGA Pipeline
# Provides structured logging with configurable levels, file output, and performance monitoring

# Load required libraries
suppressMessages({
  library(futile.logger)
  library(yaml)
})

# Global variables for logging configuration
.log_config <- NULL
.log_session_start <- NULL
.log_operation_stack <- list()

# Initialize logging system
initialize_logging <- function(config_path = "config.yaml") {
  # Load configuration
  config <- yaml.load_file(config_path)
  .log_config <<- config$logging
  
  # Set up logging directory
  logs_dir <- file.path(getwd(), config$paths$logs_dir)
  dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Configure logging level
  log_level <- switch(.log_config$level,
    "DEBUG" = DEBUG,
    "INFO" = INFO,
    "WARN" = WARN,
    "ERROR" = ERROR,
    INFO  # Default to INFO
  )
  
  # Set up console logging
  if (.log_config$console_output) {
    flog.appender(appender.console(), name = "console")
    flog.threshold(log_level, name = "console")
  }
  
  # Set up file logging
  if (.log_config$file_output) {
    log_filename <- strftime(Sys.time(), .log_config$file_pattern)
    log_filepath <- file.path(logs_dir, log_filename)
    
    flog.appender(appender.file(log_filepath), name = "file")
    flog.threshold(log_level, name = "file")
  }
  
  # Set up root logger
  flog.threshold(log_level)
  
  # Clean up old log files if rotation is enabled
  if (.log_config$rotate_logs) {
    cleanup_old_logs(logs_dir)
  }
  
  # Record session start
  .log_session_start <<- Sys.time()
  
  # Log session initialization
  log_info("=== TCGA Pipeline Logging Session Started ===")
  log_info(paste("Session ID:", format(.log_session_start, "%Y%m%d_%H%M%S")))
  log_info(paste("Log Level:", .log_config$level))
  log_info(paste("Working Directory:", getwd()))
  log_info(paste("R Version:", R.version.string))
  
  return(invisible(TRUE))
}

# Clean up old log files
cleanup_old_logs <- function(logs_dir) {
  if (!.log_config$rotate_logs) return()
  
  log_files <- list.files(logs_dir, pattern = "^tcga_pipeline_.*\\.log$", full.names = TRUE)
  
  if (length(log_files) > .log_config$max_log_files) {
    # Sort by modification time and remove oldest files
    file_info <- file.info(log_files)
    file_info <- file_info[order(file_info$mtime, decreasing = TRUE), ]
    
    files_to_remove <- rownames(file_info)[(.log_config$max_log_files + 1):nrow(file_info)]
    
    for (file_path in files_to_remove) {
      if (file.exists(file_path)) {
        file.remove(file_path)
        log_info(paste("Removed old log file:", basename(file_path)))
      }
    }
  }
}

# Core logging functions
log_debug <- function(message, ...) {
  formatted_message <- sprintf(message, ...)
  if (.log_config$console_output) flog.debug(formatted_message, name = "console")
  if (.log_config$file_output) flog.debug(formatted_message, name = "file")
}

log_info <- function(message, ...) {
  formatted_message <- sprintf(message, ...)
  if (.log_config$console_output) flog.info(formatted_message, name = "console")
  if (.log_config$file_output) flog.info(formatted_message, name = "file")
}

log_warn <- function(message, ...) {
  formatted_message <- sprintf(message, ...)
  if (.log_config$console_output) flog.warn(formatted_message, name = "console")
  if (.log_config$file_output) flog.warn(formatted_message, name = "file")
}

log_error <- function(message, ...) {
  formatted_message <- sprintf(message, ...)
  if (.log_config$console_output) flog.error(formatted_message, name = "console")
  if (.log_config$file_output) flog.error(formatted_message, name = "file")
}

# Performance monitoring functions
log_operation_start <- function(operation_name, details = NULL) {
  if (!.log_config$performance_monitoring) return(invisible(NULL))
  
  start_time <- Sys.time()
  operation_id <- paste0(operation_name, "_", format(start_time, "%H%M%S"))
  
  # Add to operation stack
  .log_operation_stack[[operation_id]] <<- list(
    name = operation_name,
    start_time = start_time,
    details = details
  )
  
  details_str <- if (!is.null(details)) paste(" -", details) else ""
  log_info("STARTED: %s%s", operation_name, details_str)
  
  return(operation_id)
}

log_operation_end <- function(operation_id, result_summary = NULL) {
  if (!.log_config$performance_monitoring || is.null(operation_id)) return(invisible(NULL))
  
  end_time <- Sys.time()
  
  if (operation_id %in% names(.log_operation_stack)) {
    operation <- .log_operation_stack[[operation_id]]
    duration <- as.numeric(difftime(end_time, operation$start_time, units = "secs"))
    
    result_str <- if (!is.null(result_summary)) paste(" -", result_summary) else ""
    log_info("COMPLETED: %s (%.2f seconds)%s", operation$name, duration, result_str)
    
    # Remove from stack
    .log_operation_stack[[operation_id]] <<- NULL
  }
  
  return(invisible(NULL))
}

# Error handling with retry logic
with_error_handling <- function(expr, operation_name = "Operation", max_retries = NULL, retry_delay = NULL) {
  # Use config defaults if not specified
  if (is.null(max_retries)) max_retries <- .log_config$error_handling$retry_attempts
  if (is.null(retry_delay)) retry_delay <- .log_config$error_handling$retry_delay
  
  operation_id <- log_operation_start(operation_name)
  
  for (attempt in 1:max_retries) {
    tryCatch({
      result <- eval(expr)
      log_operation_end(operation_id, "Success")
      return(result)
    }, error = function(e) {
      if (attempt < max_retries) {
        log_warn("%s failed (attempt %d/%d): %s. Retrying in %d seconds...", 
                 operation_name, attempt, max_retries, e$message, retry_delay)
        Sys.sleep(retry_delay)
      } else {
        log_error("%s failed after %d attempts: %s", operation_name, max_retries, e$message)
        log_operation_end(operation_id, "Failed")
        
        if (.log_config$error_handling$continue_on_error) {
          log_warn("Continuing pipeline execution despite error in %s", operation_name)
          return(NULL)
        } else {
          stop(sprintf("Pipeline stopped due to error in %s: %s", operation_name, e$message))
        }
      }
    })
  }
}

# Data quality logging
log_data_summary <- function(data_object, data_name, data_type = "matrix") {
  if (is.null(data_object)) {
    log_warn("Data object '%s' is NULL", data_name)
    return(invisible(NULL))
  }
  
  if (data_type == "matrix" || data_type == "data.frame") {
    log_info("Data Summary - %s: %d rows x %d columns", 
             data_name, nrow(data_object), ncol(data_object))
    
    if (any(is.na(data_object))) {
      na_count <- sum(is.na(data_object))
      total_count <- length(data_object)
      na_percentage <- round(na_count / total_count * 100, 2)
      log_warn("Data Quality - %s: %d missing values (%.2f%%)", 
               data_name, na_count, na_percentage)
    }
  } else if (data_type == "vector") {
    log_info("Data Summary - %s: %d elements", data_name, length(data_object))
    
    if (any(is.na(data_object))) {
      na_count <- sum(is.na(data_object))
      na_percentage <- round(na_count / length(data_object) * 100, 2)
      log_warn("Data Quality - %s: %d missing values (%.2f%%)", 
               data_name, na_count, na_percentage)
    }
  }
  
  return(invisible(NULL))
}

# Session summary logging
log_session_summary <- function() {
  if (is.null(.log_session_start)) return(invisible(NULL))
  
  session_duration <- as.numeric(difftime(Sys.time(), .log_session_start, units = "secs"))
  
  log_info("=== TCGA Pipeline Session Summary ===")
  log_info("Session Duration: %.2f seconds (%.2f minutes)", session_duration, session_duration/60)
  log_info("Session Start: %s", format(.log_session_start, "%Y-%m-%d %H:%M:%S"))
  log_info("Session End: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  
  # Log any operations still in stack (shouldn't happen in normal execution)
  if (length(.log_operation_stack) > 0) {
    log_warn("Uncompleted operations detected:")
    for (op_id in names(.log_operation_stack)) {
      op <- .log_operation_stack[[op_id]]
      log_warn("  - %s (started at %s)", op$name, format(op$start_time, "%H:%M:%S"))
    }
  }
  
  log_info("=== End of Session ===")
  return(invisible(NULL))
}

# Memory usage logging
log_memory_usage <- function(operation_name = "Memory Check") {
  if (!.log_config$performance_monitoring) return(invisible(NULL))
  
  # Get memory usage in MB
  mem_usage <- pryr::mem_used()
  mem_mb <- as.numeric(mem_usage) / 1024^2
  
  log_info("Memory Usage - %s: %.2f MB", operation_name, mem_mb)
  
  return(invisible(NULL))
}

# File operation logging
log_file_operation <- function(operation, file_path, success = TRUE, details = NULL) {
  details_str <- if (!is.null(details)) paste(" -", details) else ""
  
  if (success) {
    log_info("File %s: %s%s", toupper(operation), file_path, details_str)
  } else {
    log_error("File %s FAILED: %s%s", toupper(operation), file_path, details_str)
  }
  
  return(invisible(NULL))
}

# Progress logging for long operations
log_progress <- function(current, total, operation_name = "Processing", details = NULL) {
  if (total == 0) return(invisible(NULL))
  
  percentage <- round(current / total * 100, 1)
  details_str <- if (!is.null(details)) paste(" -", details) else ""
  
  log_info("Progress - %s: %d/%d (%.1f%%)%s", operation_name, current, total, percentage, details_str)
  
  return(invisible(NULL))
}

# Export all logging functions
export_logging_functions <- function() {
  list(
    initialize_logging = initialize_logging,
    log_debug = log_debug,
    log_info = log_info,
    log_warn = log_warn,
    log_error = log_error,
    log_operation_start = log_operation_start,
    log_operation_end = log_operation_end,
    with_error_handling = with_error_handling,
    log_data_summary = log_data_summary,
    log_session_summary = log_session_summary,
    log_memory_usage = log_memory_usage,
    log_file_operation = log_file_operation,
    log_progress = log_progress
  )
}