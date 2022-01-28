## Loads all packages and defines how to handle NAMESPACE conflicts
source("./packages.R")

## Load all R files in R/ folder
lapply(list.files("./R", full.names = TRUE), source)
pkg_list <- extract_pkg_names("packages.R")
# tar_option_set(
#   packages = c("conmat"), 
#   imports = c("conmat")
# )

# 1)
# pull out numbers of new case per day 
# detected via symptomatic surveillance 
# detected via contact tracing

# matches to data collected by jurisdictions on cases 
# and reason for test

# 2)
# check trends in proportion of infections detected via each mode of detection over time

# 3)
# testing heuristic....
# estimated number of infections = average cases found by contact tracing + average cases found by symptomatic pres * correction factor

# correction factor = 1/(symptomatic fraction * test-seeking fraction)

# case ascertainment = average cases found by contact tracing + average cases found by symptomatic pres/ estimated number of infections

# compare case ascertainment to true case ascertainment

# multiple simulations - consider TTIQ response strategies
# generate cloud of true and estimated ascertainments

# 4) 
# repeat step 3) for different specified parameter values in the ABM i.e. symptomatic fraction and test-seeking fraction

# also vary parameters when computing correction factor i.e. different to 10% relative to ABM parameters

# 5) 
# repeat steps 3) and 4) with different values of R

# and consider case where R fluctuates above and below 1 

#####################

future::plan(multisession(workers = 8))

tar_plan(
  

  
# file path example ----  
  # tar_file(
  #   cases_nsw_path, 
  #   file.path(
  #     data_path,
  #     "CASES_FROM_20200701_0000_TO_20210913_1115.xlsx"
  #   )
  # ),
  
#  cases_nsw = read_cases_nsw(cases_nsw_path),

# output file example ----
  # tar_file(scenario_test_turnaround_time_probs_path,{
  #   write_csv_return_path(scenario_test_turnaround_time_probs,
  #                         "outputs/scenario_test_turnaround_time_probs.csv")
  # }),
  
# output fig example ----
  # tar_file(plot_hist_delay_samples_v_data_path, {
  #   ggsave_write_path(
  #     plot = plot_hist_delay_samples_v_data,
  #     path = "figs/hist_delay_samples_v_data.png",
  #     width = 8,
  #     height = 8
  #   )
  # }),
  
# a doc to view plots etc ----
  # tar_render(explore, "doc/explore.Rmd", intermediates_dir="./"),
  
)
