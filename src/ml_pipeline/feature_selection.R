require(config)
require(logger)
require(tidyverse)
require(broom.helpers)
require(readr)

source("Boruta/boruta.R")

fs_boruta <- function(data_df, quantileThresh = 0.8) {
  logs <- list()
  
  set.seed(Sys.time())
  results <- Boruta(DiseaseState ~ .,
                    data = data_df,
                    pValue = 0.01,
                    quantileThresh = quantileThresh,
                    doTrace = 0)
  
  logs[['boruta_n_confirmed']] <- results$finalDecision[results$finalDecision == 'Confirmed'] %>% length
  logs[['boruta_n_rejected']] <- results$finalDecision[results$finalDecision == 'Rejected'] %>% length
  logs[['boruta_n_tentative']] <- results$finalDecision[results$finalDecision == 'Tentative'] %>% length
  
  log_debug("Boruta: {logs[['boruta_n_confirmed']]} confirmed, {logs[['boruta_n_rejected']]} rejected, {logs[['boruta_n_tentative']]} tentative.")
  
  selected_cols <- getSelectedAttributes(results, withTentative = TRUE) %>%
    .clean_backticks()  # for variable names with spaces
  
  selected_cols <- c('DiseaseState', selected_cols)
  filtered_data <- data_df %>% select(all_of(selected_cols))

  return(list(
    logs = logs,
    filtered_data = filtered_data
  ))
}


fs_utest <- function(data_df,
                     thresh_type = 'percentage',
                     thresh_pvalue = 0.1,
                     thresh_percentage = 0.2,
                     thresh_count = 30,
                     fallback = FALSE) {
  logs <- list()
  
  healthy_ds <- data_df %>% filter(DiseaseState == 'healthy')
  disease_ds <- data_df %>% filter(DiseaseState == 'disease')
  
  cols <- colnames(data_df)  
  cols <- cols[cols != 'DiseaseState']
  
  p_values <- sapply(cols, function(col_name) {
    return(wilcox.test(healthy_ds[[col_name]], disease_ds[[col_name]])$p.value)
  })
  
  results <- data.frame(col=cols, p_value=unname(p_values))
  results$corrected_p_values <- p.adjust(results$p_value, method = "fdr")
  results <- results[order(results$corrected_p_values),]
  
  if (thresh_type == 'pvalue')
    filtered_results <- results %>% filter(corrected_p_values <= thresh_pvalue)
  else if (thresh_type == 'count')
    filtered_results <- results[1:thresh_count, ]
  else if (thresh_type == 'percentage')
    filtered_results <- results[1:(as.integer(nrow(results) * thresh_percentage)), ]
  else
    stop('U test: Wrong threshold type. Available options: pvalue, count, percentage')
  
  
  logs[['utest_n_confirmed']] <- nrow(filtered_results)
  
  if (fallback) {
    min_feature_num <- as.integer(nrow(results) * thresh_percentage)
    if (nrow(filtered_results) < min_feature_num) {
      log_debug("Only {nrow(filtered_results)} features were selected. Fallback to {min_feature_num} features.")
      filtered_results <- results[1:min_feature_num, ]
    }
  }
  
  selected_cols <- c('DiseaseState', filtered_results[['col']])
  filtered_data <- data_df %>% select(all_of(selected_cols))

  return(list(
    logs = logs,
    filtered_data = filtered_data
  ))
}


fs_contribs_only <- function(data_df, curr_ds_name, curr_feature_set_type,
                             contributors_list = config::get('paths')$contributors_list) {
  logs <- list()
  
  # Read list of features previously found as contributors
  selected_cols <- read_delim(contributors_list, 
                             delim = "\t", 
                             escape_double = FALSE, 
                             trim_ws = TRUE, 
                             show_col_types = FALSE) %>%
    filter(dataset == curr_ds_name) %>%
    filter(feature_set_type == curr_feature_set_type) %>%
    pull(feature)
  
  # Keep only these features
  logs[['n_contributors_used']] <- length(selected_cols)
  selected_cols <- c('DiseaseState', selected_cols)
  filtered_data <- data_df %>% select(all_of(selected_cols))
  
  return(list(
    logs = logs,
    filtered_data = filtered_data
  ))
}

fs_select_features <- function(train_df,
                               test_df,
                               fs_type,
                               ds_name = NULL,
                               feature_set_type = NULL) {
  logs <- list()
  
  if (fs_type == 'none') {
    return(list(
      logs = logs,
      train_df = train_df,
      test_df = test_df
    ))
  } else if (fs_type == 'contribsOnly') {
    results <- fs_contribs_only(train_df, ds_name, feature_set_type)
  } else if (fs_type == 'utest') {
    results <- fs_utest(train_df)
  } else if (fs_type == 'Boruta80') {
    results <- fs_boruta(train_df, quantileThresh = 0.8)
  } else if (fs_type == 'Boruta90') {
    results <- fs_boruta(train_df, quantileThresh = 0.9)
  } else {
    stop('Wrong FS type.')
  }
  
  logs <- c(logs, results$logs)
  filtered_data <- results$filtered_data
  
  logs[['n_features_fs_selected']] <- ncol(filtered_data) - 1
  
  # If no features were selected, take all the features
  if (ncol(filtered_data) == 1) {
    filtered_data <- train_df
    log_warn("0 features selected by feature selection. Taking all features instead.")
  }
  
  # Match the test data to have the same features as the train data
  test_df <- test_df %>% select(names(filtered_data))

  return(list(
    logs = logs,
    train_df = filtered_data,
    test_df = test_df
  ))
}
