library(config)
library(tidymodels)
library(parallel)
library(optparse)

if (!exists("name__")) {
  name__ = "ml_pipeline"
}

source('utils.R')
source('preprocessing.R')
source('feature_selection.R')
source('clustering.R')
source('postprocess_results.R')

# Creates a random forest model specification (using the tidymodels framework).
# If configured to run tuning, classifier hyper-parameters are tuned here,
#  otherwise, default hyper-parameters are used.
# (Model is not trained here)
ml_create_rf_model <- function(run_name,
                               train_df,
                               should_tune) {
  logs = list()
  
  # Create a recipe
  data_recipe <- train_df %>%
    recipe(DiseaseState ~ .) %>%
    prep()
  
  set.seed(222)
  folds <- vfold_cv(train_df,
                    v = config::get('tuning_n_folds'),
                    repeats = config::get('tuning_n_repeats'))
  
  if (should_tune) {
    set.seed(333)
    
    # Create a random forest model and a workflow object
    tune_spec <- rand_forest(mtry = tune(),
                             trees = 500,
                             min_n = tune()) %>%
      set_mode("classification") %>%
      set_engine("ranger", importance = "permutation")
    
    data_workflow <- workflow() %>%
      add_recipe(data_recipe) %>%
      add_model(tune_spec)
    
    if (config::get('use_tune_grid')) {
      tune_res <- tune_grid(
        resamples = folds,
        grid = 20
      )
    } else {
      # Bayes tuning
      initial_mtry = as.integer(sqrt(ncol(train_df)))
      tune_res <- data_workflow %>%
        tune_bayes(
          resamples = folds,
          param_info = parameters(mtry(range = c(
            2, initial_mtry * 2
          )),
          min_n(range = c(1, 15))),
          iter = 30,
          metrics = metric_set(roc_auc)
        )
    }
    
    # # Graph the tuning params
    # plot <- tune_res %>%
    #   collect_metrics() %>%
    #   filter(.metric == "roc_auc") %>%
    #   select(mean, min_n, mtry) %>%
    #   pivot_longer(min_n:mtry,
    #                values_to = "value",
    #                names_to = "parameter") %>%
    #   ggplot(aes(value, mean, color = parameter)) +
    #   geom_point(show.legend = FALSE) +
    #   facet_wrap( ~ parameter, scales = "free_x") +
    #   labs(x = NULL, y = "AUC")
    # utils_save_plot(
    #   plot,
    #   sprintf('%s - Hyperparameters', run_name),
    #   width = 600,
    #   height = 300
    # )
    
    # Finalize the model and evaluate results
    best_auc <- select_best(tune_res, "roc_auc")
    logs[['mean_inner_cv_auc_best_tuning_params']] <-
      show_best(tune_res, n = 1, metric = 'roc_auc')$mean
    
    final_rf <- finalize_model(tune_spec,
                               best_auc)
  } else {
    logs[['mean_inner_cv_auc_best_tuning_params']] <- NA
    
    final_rf <- rand_forest(mode = "classification", min_n = 5) %>%
      set_mode("classification") %>%
      set_engine("ranger", importance = "permutation")
  }
  
  return(list(
    logs = logs,
    final_rf = final_rf,
    data_recipe = data_recipe
  ))
}

# Train a random forest classifier
ml_fit_model_on_train <- function(train_df, 
                                  rf_obj, 
                                  data_recipe) {
  logs = list()
  
  set.seed(555)
  final_workflow <- workflow() %>%
    add_recipe(data_recipe) %>%
    add_model(rf_obj)
  
  # Fit and document params
  final_model <- final_workflow %>% fit(data = train_df)
  fit_obj <- final_model %>% extract_fit_parsnip()
  
  logs[['rf_ntrees']] <- fit_obj$fit$num.trees
  logs[['rf_mtry']] <- fit_obj$fit$mtry
  logs[['rf_min_n']] <- fit_obj$fit$min.node.size
  
  feature_importance <- ranger::importance_pvalues(fit_obj$fit,
                                                   method = 'altmann',
                                                   formula = as.formula(data_recipe),
                                                   data = train_df)
  # Transform row names (feature name) to a separate column
  feature_importance <- feature_importance %>%
    as.data.frame() %>%
    rownames_to_column(var = "feature")
  
  return(
    list(
      logs = logs,
      fitted_wf = final_model,
      feature_importance = feature_importance
    )
  )
}


ml_evaluate_test <- function(run_name, 
                             fitted_wf, 
                             test_df) {
  logs = list()
  
  set.seed(666)
  # this will automatically bake since the recipe is in the workflow
  final_result <-
    predict(fitted_wf, new_data = test_df, type = 'prob') %>%
    bind_cols(test_df)
  
  out_of_fold_test_auc <- final_result %>%
    roc_auc(DiseaseState, .pred_disease)
  logs[['out_of_fold_test_auc']] <- out_of_fold_test_auc$.estimate
  
  plot <- final_result %>%
    roc_curve(DiseaseState, .pred_disease) %>%
    autoplot() +
    ggtitle(sprintf('%s - ROC curve', run_name))
  utils_save_plot(plot,
                  sprintf('%s - ROC curve', run_name),
                  width = 400,
                  height = 400)
  
  return(list(logs = logs))
}


ml_select_features <- function(train_df,
                               test_df,
                               fs_type,
                               ds_name = NULL,
                               feature_set_type = NULL) {
  logs <- list()
  
  # Document parameters
  logs[['h_ratio_train']] <-
    round(train_df %>% 
            filter(DiseaseState == 'healthy') %>% 
            nrow() /
            nrow(train_df), 3)
  
  logs[['h_ratio_test']] <-
    round(test_df %>% 
            filter(DiseaseState == 'healthy') %>% 
            nrow() /
            nrow(test_df), 3)
  
  log_trace("H ratio train: {logs[['h_ratio_train']]}, H ratio test: {logs[['h_ratio_test']]}")
  logs[['n_features_origin']] <- ncol(train_df) - 1
  logs[['n_features_origin_T']] <- sum(startsWith(colnames(train_df), 'T__'))
  logs[['n_features_origin_G']] <- sum(startsWith(colnames(train_df), 'G__'))
  logs[['n_features_origin_P']] <- sum(startsWith(colnames(train_df), 'P__'))
  logs[['n_features_origin_M']] <- sum(startsWith(colnames(train_df), 'M__'))
  
  # Feature selection
  results <- fs_select_features(train_df,
                                test_df,
                                fs_type,
                                ds_name,
                                feature_set_type)
  logs <- c(logs, results$logs)
  train_df <- results$train_df
  test_df <- results$test_df
  
  logs[['n_features_for_train_final']] <- ncol(train_df) - 1
  logs[['n_features_for_train_final_T']] <- sum(startsWith(colnames(train_df), 'T__'))
  logs[['n_features_for_train_final_G']] <- sum(startsWith(colnames(train_df), 'G__'))
  logs[['n_features_for_train_final_P']] <- sum(startsWith(colnames(train_df), 'P__'))
  logs[['n_features_for_train_final_M']] <- sum(startsWith(colnames(train_df), 'M__'))
  log_debug("N features original: {logs[['n_features_origin']]}, after FS: {logs[['n_features_for_train_final']]}")
  
  return(list(
    logs = logs,
    train_df = train_df,
    test_df = test_df
  ))
}


ml_ds_general_metrics <- function(data_df) {
  logs <- list()
  
  logs[['n_healthy']] <- 
    data_df %>%
    filter(DiseaseState == 'healthy') %>% 
    nrow()
  
  logs[['n_disease']] <- 
    data_df %>%
    filter(DiseaseState != 'healthy') %>% 
    nrow()
  
  logs[['h_ratio_total']] <-
    round(data_df %>% 
            filter(DiseaseState == 'healthy') %>% 
            nrow() /
            nrow(data_df), 3)
  
  log_debug(sprintf('N healthy: %d, N disease: %d',
                   logs[['n_healthy']],
                   logs[['n_disease']]))
  
  return(list(logs = logs))
}

ml_outer_fold_iteration <- function(run_name,
                                    train_df,
                                    test_df,
                                    should_tune,
                                    fs_type,
                                    ds_name = NULL) {
  logs <- list()
  log_debug('Fold iteration: {run_name}')

  # Feature selection
  results <- ml_select_features(train_df,
                                test_df,
                                fs_type,
                                ds_name,
                                feature_set_type = strsplit(run_name, " ")[[1]][5])
  
  logs <- c(logs, results$logs)
  train_df <- results$train_df
  test_df <- results$test_df
  
  # Less than 2 features (3 columns including DiseaseState) is problematic
  # for pvalue calculation
  if (ncol(train_df) < 3) {
    log_warn("Not enough features to train a model")
    return(list(logs = logs, feature_importance = list()))
  }
  
  # Choose hyperparameters and create model
  results <- ml_create_rf_model(run_name, train_df, should_tune)
  logs <- c(logs, results$logs)
  rf_obj <- results$final_rf
  data_recipe <- results$data_recipe
  
  # Fit model
  results <-
    ml_fit_model_on_train(train_df, rf_obj, data_recipe)
  logs <- c(logs, results$logs)
  fitted_wf <- results$fitted_wf
  feature_importance <- results$feature_importance
  
  # Evaluate model on test data
  results <- ml_evaluate_test(run_name, fitted_wf, test_df)
  logs <- c(logs, results$logs)
  
  return(list(logs = logs, feature_importance = feature_importance))
}

# Run the machine learning pipeline on a specific dataset
#  and a specific feature type
ml_pipeline_single_run <- function(run_name,
                                   data_df,
                                   should_cv = TRUE,
                                   should_tune = FALSE,
                                   fs_type = 'none',
                                   ds_name = NULL) {
  logs <- list()
  
  log_info('Running {run_name}')
  
  # Document general metrics on the data
  results <- ml_ds_general_metrics(data_df)
  logs <- c(logs, results$logs)
  
  
  # Run the pipeline with cross validation
  set.seed(111)
  folds <- vfold_cv(data_df,
                    v = config::get('outer_n_folds'),
                    repeats = config::get('outer_n_repeats'))
  
  # If set to run without CV take only first split
  if (!should_cv)
    folds <- folds %>% slice(1:1)
  logs[['n_folds_total']] <- nrow(folds)
  
  # Run the CV
  all_feature_importances = tibble()
  cv_results <- lapply(folds$splits, function(split) {
    fold_id <- paste(split$id, collapse = " ")
    
    # Train and evaluate a classifier on this specific split
    results <- ml_outer_fold_iteration(
      run_name = paste0(run_name, " (", fold_id, ")"),
      train_df = training(split),
      test_df = testing(split),
      should_tune = should_tune,
      fs_type = fs_type,
      ds_name = ds_name
    )
  
    return(list(
      logs = c(fold_id = fold_id, results$logs),
      feature_importance = bind_cols(results$feature_importance, fold_id = fold_id)
    ))
  })
  
  # Flat the results into a single object
  cv_results <- do.call('rbind', cv_results) %>% as_tibble
  cv_feature_importance <- bind_rows(cv_results$feature_importance)
  cv_logs <- bind_rows(cv_results$logs)
  
  # Pack the results
  results <- tibble_row(!!!logs)
  results$cv_results <- list(cv_logs)
  results$mean_out_of_fold_test_auc <- mean(cv_logs$out_of_fold_test_auc)
  results$feature_importance <- list(cv_feature_importance)

  log_info("Mean out of fold test AUC: {results$mean_out_of_fold_test_auc}")
  
  return(results)
}


ml_build_run_name <- function(params_combo) {
  
  # Extract mapping of long to short names from config file
  params_short_names <- lapply(config::get('params_short_names'), unlist)
  
  # Use mapping to get short names per pipeline setting
  run_name_l <- sapply(names(params_combo), function(p) { 
    params_short_names[[p]][as.character(params_combo[[p]])] 
  }, USE.NAMES = FALSE)
  
  # Concatenate
  run_name <- paste(run_name_l, collapse = ' ')
  
  return(run_name) 
}


ml_main <- function(ds_name) {
  log_info("Processing {ds_name}")
  
  # Create output directories
  utils_create_output_dirs()
  
  # Run the pipeline
  # Preprocess data, i.e. get feature sets
  feature_sets <- prep_preprocess_dataset(ds_name)
  
  # Get pipeline settings from config (one or more) 
  all_params_combo <- expand.grid(config::get('params_combo'))
  
  # Run pipeline
  all_results <- apply(all_params_combo, 1, function (params_combo) {
    # Get compact run name
    run_name <- ml_build_run_name(params_combo)
    
    # Set the parameters for the pipeline
    func_params = list(
      run_name = paste0(ds_name, ' - ', run_name),
      data_df = feature_sets[[params_combo[['feature_set_type']]]],
      should_tune = params_combo[['should_tune']],
      fs_type = params_combo[['fs_type']],
      should_cv = TRUE,
      ds_name = ds_name
    )
    
    results <- do.call('ml_pipeline_single_run', func_params)
    
    results$run_name = run_name
    results$feature_set_type = params_combo[['feature_set_type']]
    results$should_tune = params_combo[['should_tune']]
    results$fs_type = params_combo[['fs_type']]
    return(results)
  })
  all_results <- bind_rows(all_results)
  
  # Save a CSV file with the folds results
  result_path <- sprintf(config::get('paths_templates')$cv_results_csv, ds_name)
  
  all_results %>% 
    mutate(dataset = ds_name) %>%
    relocate(dataset) %>%
    select(-feature_importance) %>%
    unnest(cols = cv_results) %>%
    mutate(run_name = paste0('+"', run_name, '"')) %>%  # so that excel won't show this as a formula
    utils_save_tsv_table(result_path, sep = ',')
  
  # Save a CSV file with features importance
  feature_importance_path <- sprintf(config::get('paths_templates')$feature_importance_csv, ds_name)
  
  all_results %>% 
    mutate(dataset = ds_name) %>%
    relocate(dataset) %>%
    select(c(dataset, run_name, feature_importance)) %>%
    unnest(cols = feature_importance) %>%
    mutate(run_name = paste0('+"', run_name, '"')) %>%  # so that excel won't show this as a formula
    utils_save_tsv_table(`feature_importance_path`, sep = ',')
}

# Parse command line arguments
ml_parse_args <- function() {
  option_list <- list(
    make_option(c("-d", "--dataset"), type = "character", default = NULL, 
                help = "dataset name", metavar = "character")
  ); 
  
  opt_parser <- OptionParser(option_list = option_list);
  opts <- parse_args(opt_parser);
  return(opts)
}

###################################
# Main 
###################################

if (name__ == "ml_pipeline") {
  # Set working directory 
  setwd(getwd())
  
  # Set log level (globally)
  log_threshold(config::get('log_level'))
  
  # Parse arguments (will take default values if none given)
  args <- ml_parse_args()
  
  # The ML pipeline could be executed in a few different ways:
  # Option 1: User specifies a dataset to run on using the "-d" command line argument
  # ----------------------------------------------------------------------------------
  if (!is.null(args$dataset)) {
    
    ml_main(args$dataset)
    
    # Option 2: Run in test mode (i.e. use 'test' configuration)
    # ----------------------------------------------------------------------------------
  } else {
    
    Sys.setenv(R_CONFIG_ACTIVE = "test")  # Set to test configuration mode globally (i.e. use 'test' configurations in config.yml)
    log_info('Running with test configuration...')
    ml_main('crc_wang') 
    ml_main('crc_zeller') 
    ml_main('t1d_alkanani') 
    ml_main('crc_zeller_2014') 
    ml_main('adenomas_feng_2015') 
    ml_main('esrd_wang_2020') 
    ml_main('uc_franzosa_2019')
    ml_main('cd_franzosa_2019')
    post_prepare_rdata()
    post_save_contributors()
  }
  
  # Option 3: Run on a list of datasets/settings in parallel. 
  # For this use the 'run_pipeline_parallel.R' script.
}
