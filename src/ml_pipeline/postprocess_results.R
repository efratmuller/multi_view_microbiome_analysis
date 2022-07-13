require(config)
require(metap)
require(readr)
require(dplyr)

source('utils.R')

# Combine features' importance p values using Fisher's method (function for use in following command)
post_combine_p_vals_fisher <- function(pvals) {
  if (length(pvals) == 1) return(pvals)
  return(sumlog(pvals)$p)
}

# Combine all cv_results (saved as one file per dataset) into a single table, 
#  re-format some columns and add extra information.
# Each row in the output represents a single fold/repeat in a specific dataset 
#  using a specific setting ("run_name").
# Also writes the table to a file by default.
post_load_cv_results <- function(dir_path = config::get("paths")$results_tables_dir, 
                                 output_file = config::get("paths")$combined_cv_results) {
  files <- list.files(dir_path, pattern = "\\_pipeline.csv$")
  all_data <- bind_rows(lapply(files, function(result_file_name) {
    return(read.csv(file.path(dir_path, result_file_name)))
  }))
  
  # Remove the extra plus at the beginning, added by ml_pipeline.R for excel
  all_data$run_name <- substring(all_data$run_name, 2, nchar(all_data$run_name))
  
  # Order categories for consistent plotting
  all_data <- all_data %>%
    mutate(run_name = factor(run_name)) %>%
    mutate(feature_set_type = 
             factor(feature_set_type, 
                    levels = names(config::get("params_short_names")$feature_set_type)))
  
  # Mark the shotgun studies using curatedMetagenomicData list of included studies
  shotgun_studies <- read_csv(config::get("paths")$curatedMetagenomicData_study_list, 
                              comment = "#", 
                              show_col_types = FALSE)
  all_data$metagenomics_type <- ifelse(all_data$dataset %in% 
                                         c(shotgun_studies$new_study_name, config::get("shotgun_datasets_not_in_CMD")), 
                                       "Shotgun", 
                                       "16S")
  
  # Write to file
  if (! is.null(output_file))
    utils_save_tsv_table(all_data, output_file, sep = ',')
  
  return(all_data)
}


# Combine all feature importance results (saved as one file per dataset) into a single table, 
#  re-format some columns and add extra information.
# Also writes the table to a file by default.
post_load_feature_importance <- function(dir_path = config::get("paths")$results_tables_dir, 
                                         output_file = config::get("paths")$combined_feature_importance) {
  files <- list.files(dir_path, pattern = "\\_feature_importance.csv$")
  all_data <- bind_rows(lapply(files, function(result_file_name) {
    return(read.csv(file.path(dir_path, result_file_name)))
  }))
  
  # Remove the extra plus at the beginning, added by ml_pipeline.R for excel
  all_data$run_name <- substring(all_data$run_name, 2, nchar(all_data$run_name))
  
  # Join extra data from cv results
  cv_results <- post_load_cv_results(dir_path, output_file = NULL) %>%
    select(c(dataset, run_name, mean_out_of_fold_test_auc)) %>%
    distinct()
  all_data <- all_data %>%
    inner_join(cv_results, by = c('dataset', 'run_name'))
  
  # Add flag if feature was part of a cluster
  all_data <- all_data %>% 
    mutate(is_cluster_rep = grepl("_Cluster", feature)) %>%
    mutate(cluster_id = ifelse(is_cluster_rep, gsub("^.*_Cluster", "", feature), NA))
  
  # Add feature type
  all_data <- all_data %>% 
    mutate(feature_type = substr(feature, 1, 1))
  
  # Add per feature the number of times it was selected and other summary stats
  #  (min = 1, max = 50 (5 folds 10 repeats). Features never selected are not in the list)
  feat_imp_sum <- all_data %>% 
    group_by(dataset, run_name, feature, feature_type, is_cluster_rep, cluster_id) %>% 
    summarise(n_times_selected = n(), 
              mean_importance = mean(importance),
              combined_p = post_combine_p_vals_fisher(pvalue),
              .groups = "drop") %>%
    ## Add FDR
    group_by(dataset, run_name) %>% 
    mutate(combined_fdr = p.adjust(combined_p, method = "fdr")) %>%
    ungroup()
  
  # Add these summary stats to main table
  all_data <- all_data %>%
    left_join(feat_imp_sum, 
              by = c("dataset",
                     "run_name",
                     "feature",
                     "feature_type",
                     "is_cluster_rep",
                     "cluster_id"))
  
  # Map compound id's to names
  kegg_metab_names_map <- utils_get_kegg_compound_names_mapping()
  all_data <- all_data %>%
    mutate(mtb_only = ifelse(grepl("^M__", feature), 
                             gsub("M__","",gsub("_Cluster[0-9]*$","",feature)),
                             NA)) %>%
    mutate(pretty_metab_name = unname(kegg_metab_names_map[mtb_only])) %>%
    # After mapping, paste back the 'M_' prefix and cluster suffix where relevant
    mutate(pretty_metab_name = ifelse(!is.na(pretty_metab_name),
                                      paste0("M__", pretty_metab_name),
                                      NA)) %>%
    mutate(pretty_metab_name = ifelse(!is.na(cluster_id) & !is.na(pretty_metab_name),
                                      paste0(pretty_metab_name, "_Cluster", cluster_id),
                                      pretty_metab_name)) %>%
    mutate(pretty_feature_name = coalesce(pretty_metab_name, feature)) %>%
    select(-pretty_metab_name, -mtb_only) %>%
    # Set some maximal string length for plotting purposes
    mutate(pretty_feature_name = ifelse(nchar(pretty_feature_name) > 40,
                                        paste0(substr(pretty_feature_name, 1, 40),"..."),
                                        pretty_feature_name))
  
  # TODO: also map MetaCyc pathway names, from: config::get('paths')$metacyc_pathways_names
  
  # Write to file
  if (! is.null(output_file))
    utils_save_tsv_table(all_data, output_file)
  
  return(all_data)
}

post_summarize_feat_imp <- function(feat_imp,
                                    fdr_threshold = 0.05,
                                    n_times_selected_threshold = 0.5,
                                    n_total_folds = config::get("outer_n_folds") * config::get("outer_n_repeats")) {
  feat_imp_sum <- feat_imp %>%
    select(dataset, run_name, 
           mean_out_of_fold_test_auc, 
           feature, feature_type, 
           mean_importance, combined_p, 
           combined_fdr, n_times_selected, pretty_feature_name, 
           is_cluster_rep, cluster_id) %>%
    distinct() %>%
    mutate(contributor = (combined_fdr <= fdr_threshold) & 
             (n_times_selected/n_total_folds > n_times_selected_threshold))
  return(feat_imp_sum)
}

# Save contributing features to a file 
post_save_contributors <- function(output_file = config::get("paths")$contributors_list) {
  # Load feature importance
  feat_imp <- post_load_feature_importance(output_file = NULL)
  
  # Summarize feature importance
  feat_imp_sum <- post_summarize_feat_imp(feat_imp) %>%
    mutate(feature_set_type = gsub("^.* ", "", run_name))
  
  # Save thin list to file
  contrib_list <- feat_imp_sum %>%
    filter(contributor) %>%
    select(dataset, feature_set_type, feature) %>%
    utils_save_tsv_table(output_file)
}


post_add_microbiomeHD_results <- function(results) {
  microbiomeHD_auc = tibble(dataset = c("cdi_schubert",
                                        "noncdi_schubert",
                                        "crc_baxter",
                                        "edd_singh",
                                        "ob_goodrich",
                                        "ob_turnbaugh",
                                        "ob_zupancic",
                                        "t1d_alkanani",
                                        "par_scheperjans",
                                        "asd_son",
                                        "crc_wang",
                                        "crc_zeller"),
                            microbiomeHD_auc = c(0.99,
                                                 0.98,
                                                 0.77,
                                                 0.96,
                                                 0.67,
                                                 0.84,
                                                 0.44,
                                                 0.71,
                                                 0.67,
                                                 0.39,
                                                 0.9,
                                                 0.82)
  )
  joined_data <- results %>%
    left_join(microbiomeHD_auc, by = 'dataset')
  return(joined_data)
}

# Prepares an rdata object with all results in it, as well as saves text files with combined results (combined over all datasets)
post_prepare_rdata <- function(dir_path = config::get("paths")$results_tables_dir, 
                               output_file = config::get("paths")$results_rdata, 
                               auc_threshold = 0.7) {
  # 1. Load ML modelling results
  ######################################
  
  ## Load (and combine into one table) all CV results 
  cv_results <- post_load_cv_results()
  
  ## Rounded AUC for plotting only
  cv_results <- cv_results %>%
    mutate(oof_test_auc_rounded = ifelse(out_of_fold_test_auc < 0.5, 0.5, out_of_fold_test_auc))
  
  # 2. Load feature importance results 
  ######################################
  
  feat_imp <- post_load_feature_importance()
  feat_imp_sum <- post_summarize_feat_imp(feat_imp)
  
  # Add feature_set_type column
  feat_imp$feature_set_type <- gsub("^.* ", "", feat_imp$run_name)
  feat_imp_sum$feature_set_type <- gsub("^.* ", "", feat_imp_sum$run_name)
  
  # 3. Remove low-performing datasets
  ######################################
  
  datasets_to_analyze <- cv_results %>%
    filter(feature_set_type == "T+G+P+M") %>%
    filter(fs_type %in% c("Boruta90","contribsOnly")) %>%
    group_by(dataset) %>%
    filter(max(mean_out_of_fold_test_auc) >= auc_threshold) %>%
    pull(dataset) %>%
    unique()
  
  datasets_discarded <- setdiff(unique(cv_results$dataset), 
                                datasets_to_analyze)
  
  message("Only the following datasets will be further analyzed:",
          paste(datasets_to_analyze, collapse = ", "))
  message("Discarded datasets:",
          paste(datasets_discarded, collapse = ", "))
  
  cv_results <- cv_results %>%
    filter(dataset %in% datasets_to_analyze)
  
  feat_imp <- feat_imp %>%
    filter(dataset %in% datasets_to_analyze) 
  
  feat_imp_sum <- feat_imp_sum %>%
    filter(dataset %in% datasets_to_analyze)

  # 4. Load clusters data
  ######################################
  
  all_clusters <- bind_rows(
    lapply(list.files(dir_path, pattern = "_clusters.csv$"),
           function (result_file_name) {
             df <- read.csv(file.path(dir_path, result_file_name))
             df$dataset <- str_split(result_file_name, "_clusters.csv")[[1]][1]
             return(df)
           })) %>%
    filter(dataset %in% datasets_to_analyze) %>%
    rename(feature_orig = feature) %>%
    rename(feature = representative)
  
  # 5. Save
  ######################################
  
  save(cv_results, 
       feat_imp, 
       feat_imp_sum, 
       all_clusters, 
       datasets_to_analyze, 
       datasets_discarded,
       file = output_file)
}



# post_ds_plot_final_results <- function(row, ds_name) {
#   ds_results <- row$data[[1]]
#   ds_name <- ds_name$dataset
#   
#   # Nest what's necessary
#   keep_cols = c('run_name', 'mean_out_of_fold_test_auc',
#                 'feature_set_type', 'should_tune', 'fs_type')
#   nest_cols = colnames(row$data[[1]])
#   nest_cols = nest_cols[!(nest_cols %in% keep_cols)]
#   
#   ds_results <- ds_results %>% nest(data = nest_cols)
#   
#               
#   auc_data <- tibble(desc = character(),
#                      test_auc = numeric(),
#                      should_tune = numeric(),
#                      fs_type = character(),
#                      feature_set_type = character(),
#                      label = character())
#   
#   auc_data <- auc_data %>% add_row(desc = 'MicrobiomeHD',
#                                    test_auc = row$microbiomeHD_auc,
#                                    should_tune = 0,
#                                    fs_type = 'none',
#                                    feature_set_type = 'MicrobiomeHD',
#                                    label = '')
#   auc_data$test_auc[is.na(auc_data$test_auc)] <- 0  # in case it's not MicrobiomeHD dataset
#   
#   for (j in 1:nrow(ds_results)) {
#     results_row <- ds_results %>%
#       slice(j:j)
#     
#     num_all_features <- mean(results_row$data[[1]]$n_features_post_clustering)
#     if (results_row$fs_type == 'none') {
#       num_features <- mean(results_row$data[[1]]$n_features_post_clustering)
#     } else {
#       num_features <- mean(results_row$data[[1]]$n_features_fs_selected)
#     }
#     
#     label = paste0(round(num_features), '/', round(num_all_features))
#     
#     auc_data <- auc_data %>% add_row(desc = results_row$run_name,
#                                      test_auc = results_row$mean_out_of_fold_test_auc,
#                                      should_tune = results_row$should_tune,
#                                      fs_type = results_row$fs_type,
#                                      feature_set_type = results_row$feature_set_type,
#                                      label = label)
#     
#   }
#   auc_data$desc <- factor(auc_data$desc, levels = auc_data$desc)  # for order
#   auc_data$should_tune <- factor(auc_data$should_tune)
#   
#   # Calculate the dash line groups
#   num_in_each_group <- auc_data %>%
#     filter(feature_set_type != 'MicrobiomeHD') %>%
#     group_by(feature_set_type) %>%
#     count
#   num_groups <- nrow(num_in_each_group)
#   num_in_each_group <- num_in_each_group$n[[1]]
#   lines_seq <- seq(from = 1.5,
#                    to = (1.5 + (num_in_each_group * (num_groups - 1))),
#                    by = num_in_each_group)
#   
#   plot_title <- sprintf('%s - Results summary', ds_name)
#   plot <- ggplot(data = auc_data, aes(colour = fs_type)) +
#     geom_point(aes(x = desc, y = test_auc), shape = 16, size = 3) +
#     geom_label(aes(x = desc, y = test_auc, label = label), nudge_y = 0.01, hjust = 'left') +
#     #geom_point(aes(x = desc, y = folds_auc), shape = 21, size = 3) +
#     scale_color_manual(values = c('none' = '#878787',
#                                   'utest' = '#00ccff',
#                                   'Boruta90' = '#ff8000',
#                                   'Boruta80' = '#b3a100')) +
#     geom_vline(xintercept = lines_seq, 
#                linetype = "dashed", colour = "gray", size = 0.5) +
#     ylim(0.5, 1) +
#     ylab('AUC') + xlab('') +
#     coord_flip() +
#     ggtitle(plot_title) + 
#     theme_bw() +
#     theme(legend.position = 'none',
#           plot.title = element_text(face="bold", size = rel(1.2), hjust = 0.5),
#           text=element_text(size = 12))
#   #print(plot)
#   utils_save_plot(plot, plot_title, width = 400, height = 1300)
#   
#   
#   # Plot summary for ds type
#   group_max_auc <- auc_data %>%
#     group_by(feature_set_type) %>%
#     summarise(max_auc = max(test_auc))
#   group_max_auc$feature_set_type <- factor(group_max_auc$feature_set_type, levels=c('MicrobiomeHD', 'T', 'P', 'M', 'T+P', 'T+M', 'T+P+M'))  # for order
#   
#   plot_title <- sprintf('%s - Best results', ds_name)
#   plot <- ggplot(data = group_max_auc, aes(x = feature_set_type, y = max_auc, colour = feature_set_type)) +
#     geom_point(size = 3) +
#     ylim(0.5, 1) +
#     ylab('AUC') + xlab('') +
#     coord_flip() +
#     ggtitle(plot_title) + 
#     theme_bw() +
#     theme(legend.position = 'none',
#           plot.title = element_text(face="bold", size = rel(1.2), hjust = 0.5),
#           text = element_text(size = 12))
#   #print(plot)
#   utils_save_plot(plot, plot_title, width = 350, height = 250)
# }


# post_summarize_best_config <- function(results) {
#   auc_per_config <- results %>%
#     unnest(cols = c(data)) %>%
#     select(c(dataset,
#              mean_out_of_fold_test_auc,
#              feature_set_type,
#              should_tune,
#              fs_type,
#              cluster_type)) %>%
#     distinct()  # discard the folds values
#     
#   best_config <- auc_per_config %>%
#     group_by(dataset, feature_set_type, cluster_type) %>%
#     filter(mean_out_of_fold_test_auc == max(mean_out_of_fold_test_auc)) %>%
#     slice(1) %>%  # discard multiple values with the same best AUC
#     select(-c(mean_out_of_fold_test_auc))
#   
#   utils_save_tsv_table(best_config, config::get('paths')$best_config_output, sep = ',')
# }


# post_plot_feature_importance <- function(feature_importance_data, group_data) {
#   # Take the most important features, they are already sorted by let's do it anyway
#   sliced_features <- feature_importance_data %>%
#     arrange(-mean_importance) %>%
#     slice(1:30)
#   
#   # Check if there was resolving
#   if (!sliced_features %in% colnames(sliced_features)) {
#     sliced_features <- sliced_features %>%
#       mutate(pretty_name = feature) %>%
#       mutate(feature_type = 'N')
#   }
#   
#   # Trim names
#   sliced_features$pretty_name <- sapply(sliced_features$pretty_name, function(feature_name) {
#     max_feature_name_length = 20
#     if (nchar(feature_name) <= max_feature_name_length) {
#       return(feature_name)
#     }
#     return(paste0(substring(feature_name, 1, max_feature_name_length), '...'))
#   })
#   
#   # Remove duplicates and factor to keep the original order
#   sliced_features <- sliced_features %>%
#     distinct(across(pretty_name), .keep_all = TRUE)
#   sliced_features$pretty_name <- factor(sliced_features$pretty_name, levels = sliced_features$pretty_name) %>%
#     fct_rev()
#   sliced_features$feature_type <- factor(sliced_features$feature_type)  # for color
#   
#   plot_title <- 
#   plot <- ggplot(data = sliced_features) +
#     geom_bar(aes(x = pretty_name, y = mean_importance, fill = feature_type),
#              stat = 'identity', show.legend = FALSE) +
#     scale_fill_manual(values = c('T' = '#66dd6e', 'P' = '#f8766d', 'M' = '#00bfc4')) +
#     coord_flip() +
#     xlab('Feature') + 
#     ylab('Importance') +
#     ggtitle(sprintf('%s - %s (AUC %1.2f)',
#                     group_data$dataset,
#                     group_data$run_name,
#                     group_data$mean_out_of_fold_test_auc)) +
#     theme_bw() 
#     theme(legend.position = "none",
#           plot.title = element_text(face="bold", size = rel(1.2), hjust = 0.5))
#   
#   print(plot)
#   utils_save_plot(plot,
#                   sprintf('%s - %s - Feature importance', group_data$dataset, group_data$run_name),
#                   width = 500, height = 600)
# }


# post_main <- function(cv_plots = FALSE,
#                       best_config = TRUE,
#                       feature_importance = FALSE) {
#   # Plot graphs for each dataset
#   if (cv_plots) {
#     cv_results <- post_load_cv_results(config::get('paths')$results_tables_dir)
#     cv_results_with_microbiomeHD <- post_add_microbiomeHD_results(cv_results)
#     #cv_results_with_microbiomeHD <- cv_results_with_microbiomeHD %>% filter(dataset == 'crc_baxter')
#     
#     group_walk(cv_results_with_microbiomeHD, post_ds_plot_final_results)
#   }
#   
#   # Create a CSV with a summary for the best config for each dataset type
#   if (best_config) {
#     cv_results <- post_load_cv_results(config::get('paths')$results_tables_dir)
#     post_summarize_best_config(cv_results)
#   }
#   
#   # Plot graphs for feature importance
#   if (feature_importance) {
#     feature_importance <- post_load_feature_importance(config::get('paths')$results_tables_dir)
#     feature_importance %>%
#       group_by(dataset, run_name, mean_out_of_fold_test_auc) %>%
#       group_walk(post_plot_feature_importance)
#   }
# }

###################################
# Main 
###################################

# dir_path <- config::get('paths')$results_tables_dir
# dummy <- post_load_cv_results(dir_path)
# dummy <- post_load_feature_importance(dir_path)


