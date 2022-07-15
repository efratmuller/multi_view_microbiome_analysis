require(config)
require(logger)
require(dplyr)

source('pathways_db.R')


prep_add_col_prefix <- function(df, prefix) {
  sample_id_col <- df %>% select(sample_id__)
  df <- df %>% select(-sample_id__)
  colnames(df) <- paste0(prefix, "__", colnames(df))
  df <- bind_cols(sample_id_col, df)
  return(df)
}


prep_load_metadata <- function(path) {
  valid_labels = c('H', config::get('disease_labels'))

    metadata_df <- read.csv(file = path, sep = '\t', header = TRUE) %>%
      as_tibble() %>%
      rename(sample_id__ = 1) %>%
      mutate(sample_id__ = as.character(sample_id__)) %>%
      select(sample_id__, DiseaseState) %>%
      filter(DiseaseState %in% valid_labels) %>%
      mutate(DiseaseState = replace(DiseaseState, DiseaseState == 'H', 'healthy')) %>%
      mutate(DiseaseState = replace(DiseaseState, DiseaseState != 'healthy', 'disease')) %>%
      mutate_at(vars(DiseaseState), factor)
  return(metadata_df)
}


prep_load_pathways <- function(path) {
  # Load pathways and filter redundant (superpathways) pathways
  pathways_df <- read.table(file = path, sep = '\t',
                            header = TRUE,
                            quote = "",
                            check.names = FALSE) %>%
    filter_pathways_df()
  
  # Transpose pathways DF and correct column names to pathways names
  pathways_df <- pathways_df %>%
    t() %>%
    data.frame() %>%
    tibble::rownames_to_column() %>%
    as_tibble()

  names(pathways_df) <- unlist(pathways_df[1, ])
  
  pathways_df <- pathways_df[-(1:2), ]
  pathways_df <- rename(pathways_df, sample_id__ = 1)
  pathways_df <- prep_add_col_prefix(pathways_df, 'P')
  
  pathways_df <- prep_convert_to_relative_abundances(pathways_df)
  return(pathways_df)
}

prep_load_metagenome <- function(path) {
  # Read table, remove "description" column if exists (outputted by picrust)
  df <- read.csv(file = path, header = TRUE, check.names = FALSE, sep = '\t') %>%
    select(-any_of(c('description'))) 
  
  if ("function" %in% names(df)) 
    df <- df %>% rename(`KO` = `function`)
  
  df <- df %>%
    tibble::column_to_rownames(var = "KO") %>%
    t() %>%
    data.frame() %>%
    tibble::rownames_to_column(var = "sample_id__") %>%
    prep_add_col_prefix('G')
  
  df <- prep_convert_to_relative_abundances(df)
  return(df)
}

prep_load_taxonomy <- function(path) {
  df <- read.csv(file = path, sep = '\t', header = TRUE, check.names = FALSE) %>%
    select(-any_of(config::get('unknown_taxa'))) %>%
    rename(sample_id__ = 1) %>%
    prep_add_col_prefix('T')
  
  df <- prep_convert_to_relative_abundances(df)
  return(df)
}


prep_load_metabolites <- function(path) {
  metabolites_df <- read.csv(file = path, sep = '\t', header = TRUE, check.names = FALSE) %>%
    rename(sample_id__ = 1) %>%
    mutate(sample_id__ = as.character(sample_id__)) %>%  # TODO: remove?
    prep_add_col_prefix('M')
  return(metabolites_df)
}


prep_join_metadata <- function(dataset_df, metadata_df) {
  all_data <- dataset_df %>%
    inner_join(metadata_df, by = 'sample_id__')
  return(all_data)
}


prep_join_features <- function(df1, df2) {
  joined_df <- df1 %>% 
    inner_join(df2, by = 'sample_id__') %>%
    relocate(sample_id__)
  return(joined_df)
}


prep_convert_to_relative_abundances <- function(df) {
  sample_id_col <- df %>%
    select(sample_id__) %>%
    mutate(sample_id__ = as.character(sample_id__))
  
  df <- df %>%
    select(-sample_id__) %>%
    mutate_all(as.numeric) %>%
    mutate(totalSum = rowSums(across(everything()))) %>%
    mutate(across(.cols = -c(totalSum), ~ .x / totalSum)) %>%
    select(-totalSum) %>%
    bind_cols(sample_id_col)
  return(df)
}


prep_sanitize_dataset <- function(df, feature_set, rare_feature_cutoff = 0.2) {
  original_nfeatures <- ncol(df)
  
  # Remove constant features
  df <- df %>% select(where(~ n_distinct(.) > 1))
  
  # Remove rare features (have <20% non-zero values)
  # Note: this means we use a prevalence filter only, regardless of abundance
  non_zero_percentage <- colSums(df != 0) / nrow(df)
  rare_features <- names(non_zero_percentage[non_zero_percentage <= rare_feature_cutoff])
  df <- df %>% select(-all_of(rare_features))
  
  log_debug(sprintf('Sanitizer: removed %d/%d (%d%%) rare/constant features for dataset %s',
                    (original_nfeatures - ncol(df)),
                    original_nfeatures,
                    as.integer(100 * (original_nfeatures - ncol(df)) / original_nfeatures),
                    feature_set))
  return(df)
}


# There is an option for override configuration clustering, used by cross.
# If NULL is supplied, the clustering config is used. Otherwise, the clustering
# specifid by "override_clustering" arg is used.
prep_preprocess_dataset <- function(ds_name, cluster_type = config::get('cluster_type')) {
  # Get file paths
  metadata_path <- sprintf(config::get('paths_templates')$metadata, ds_name)
  taxonomy_path <- sprintf(config::get('paths_templates')$taxonomy, ds_name)
  pathways_path <- sprintf(config::get('paths_templates')$pathways, ds_name)
  metabolites_path <- sprintf(config::get('paths_templates')$metabolites, ds_name)
  genes_path <- sprintf(config::get('paths_templates')$metagenome, ds_name)
  clusters_output <- sprintf(config::get('paths_templates')$clusters, ds_name)
  
  # Initialize clusters file
  cluster_init_clusters_table(clusters_output)

  # Load datasets
  metadata_df <- prep_load_metadata(metadata_path)
  taxonomy_df <- prep_load_taxonomy(taxonomy_path)
  pathways_df <- prep_load_pathways(pathways_path)
  metabolites_df <- prep_load_metabolites(metabolites_path)
  genes_df <- prep_load_metagenome(genes_path)
  
  # Run processing steps that are common to all datasets and independent of dataset combinations
  datasets = list(`T` = taxonomy_df,
                  `G` = genes_df,
                  `P` = pathways_df,
                  `M` = metabolites_df)
  
  clean_dfs <- lapply(seq_along(datasets), function(i) {
    dataset <- datasets[[i]]
    dataset_name <- names(datasets)[i]

    # 1. Sort by sample ID so the folds will be the same across different feature types
    dataset <- dataset[order(dataset$`sample_id__`),]
    
    # 2. Remove rare/constant features
    dataset <- prep_sanitize_dataset(dataset, dataset_name)
    
    # 3. Fix column names
    # (Some column names may be problematic for formulas, we thus fix them)
    colnames(dataset) <- make.names(colnames(dataset))
    return(dataset)
  })
  names(clean_dfs) <- names(datasets)
  
  # Create combinations
  t_and_g <- prep_join_features(clean_dfs$`T`, clean_dfs$`G`)
  t_and_p <- prep_join_features(clean_dfs$`T`, clean_dfs$`P`)
  t_and_m <- prep_join_features(clean_dfs$`T`, clean_dfs$`M`)
  t_g_p   <- prep_join_features(t_and_g, clean_dfs$`P`)
  t_g_m   <- prep_join_features(t_and_g, clean_dfs$`M`)
  t_p_m   <- prep_join_features(t_and_p, clean_dfs$`M`)
  t_g_p_m <- prep_join_features(t_g_p, clean_dfs$`M`)
  
  # Organize all combinations in a single list (override previous list)
  datasets = list(`T` = clean_dfs$`T`,
                  `G` = clean_dfs$`G`,
                  `P` = clean_dfs$`P`,
                  `M` = clean_dfs$`M`,
                  `T+G` = t_and_g,
                  `T+P` = t_and_p,
                  `T+M` = t_and_m,
                  `T+G+P` = t_g_p,   
                  `T+G+M` = t_g_m,   
                  `T+P+M` = t_p_m,   
                  `T+G+P+M` = t_g_p_m)
  
  # Run these final processing steps for each data combination
  finalized_dfs <- lapply(seq_along(datasets), function(i) {
    dataset <- datasets[[i]]
    dataset_name <- names(datasets)[i]
    
    # 1. Add disease state (this will be the label we wish to predict)
    dataset <- prep_join_metadata(dataset, metadata_df)
    dataset <- dataset %>% select(-`sample_id__`)
    
    # 2. Cluster highly-redundant features 
    # (We randomly choose a single representative, and document identified clusters)
    dataset <- cluster_cluster_features(
      dataset = dataset, 
      feature_set_type = dataset_name,
      cluster_type = cluster_type,
      clusters_output = clusters_output
    )
    
    return(dataset)
  })
  names(finalized_dfs) <- names(datasets)
  
  return(finalized_dfs)
}
