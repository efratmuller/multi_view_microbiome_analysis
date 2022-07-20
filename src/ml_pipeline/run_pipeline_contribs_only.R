########################################################################
# This script is meant to check how well our models perform when 
#  trained only on the list of "contributing features".
# 
# We define "contributing features" as features selected by the feature
#  selection process at least 50% of the time (over repeated folds),
#  and having an importance p value below a threshold.
# 
# This script can be run only after ml_pipeline ran on all datasets and 
#  functions from postprocess_results.R were used to combine all results 
#  and save a list of contributors only (per dataset, per run_name.. 
#  saved in config::get("paths")$contributors_list).
########################################################################

# Patch to load ml_pipeline functions without running the "main"
if (!exists("name__")) {
  name__ = "ml_pipeline_contrib_only"
}
source('ml_pipeline.R')

# Run pipeline using contributors only for all available datasets
# TODO: parallelize
# TODO: no need to re-write the clusters file (for now marked with REDUNDANT suffix)
run_pipeline_contribs_only <- function(dataset_names = NULL) {
  if (is.null(dataset_names)) 
    dataset_names <- utils_get_available_ds()
  dummy <- lapply(dataset_names, ml_main)
}

###################################
# Main 
###################################

# Set working directory 
setwd(getwd())

# Set log level (globally)
log_threshold(config::get('log_level'))

# Set "contributors" mode
Sys.setenv(R_CONFIG_ACTIVE = "contributors")

# Run pipeline on all datasets, by default - all datasets found in the ml_input package
if (name__ == "ml_pipeline_contrib_only") {
  run_pipeline_contribs_only()
}

# For testing
if (FALSE) {
  run_pipeline_contribs_only(
    c('crc_zeller', 
      't1d_alkanani', 
      'crc_zeller_2014', 
      'adenomas_feng_2015', 
      'esrd_wang_2020', 
      'uc_franzosa_2019',
      'cd_franzosa_2019')
  )
  post_prepare_rdata()
}

