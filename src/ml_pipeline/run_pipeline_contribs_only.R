########################################################################
# This script is meant to check how well our models performed when 
#  trained only on the list of contributing features found by the full
#  pipeline.
# This script can be run only after ml_pipeline ran on all datasets and 
#  functions from postprocess_results.R were used to combine all results 
#  and save a list of contributors only (per dataset, per run_name.. 
#  saved in config::get("paths")$contributors_list).
# If these models perform comparably good as the previous ones, we infer
#  that our list of contributors indeed covers all informative features.
########################################################################

# Patch to load ml_pipeline functions without running the main function there
if (!exists("name__")) {
  name__ = "ml_pipeline_contrib_only"
}
source('ml_pipeline.R')

# Run pipeline using contributors only for all datasets available
# TODO: parallelize
# TODO: no need to re-write the clusters file (for now marked with REDUNDANT suffix)
run_pipeline_contribsOnly_in_parallel <- function() {
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
  run_pipeline_contribsOnly_in_parallel()
}

if (FALSE) {
  # For testing
  ml_main('t1d_alkanani')
  ml_main('crc_wang') 
  ml_main('crc_zeller') 
  ml_main('crc_zeller_2014') 
  ml_main('adenomas_feng_2015') 
  ml_main('esrd_wang_2020') 
  ml_main('uc_franzosa_2019')
  ml_main('cd_franzosa_2019')
  post_prepare_rdata()
}

