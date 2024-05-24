# Script to generate the first_lastgene.csv
# read & unzip _nat.csv.gzip, extract the data of the first and last generation, and save it a csv file.

library(tidyverse)
library(data.table)

#Rep_n <- 10
G_n <- 5000
#P_n <- 10000
#T_n <- 500

# set the suffix of the target directories as a vector e.g., "_20220322151600", "_20220322151601", ...
dir_pattern <- c("_20220322151600") # an example

# The function which extracts the data of the first and last generation
extract_FirstLast <- function(i_file) {
  
  print(paste0("The target file is ", i_file))
  
  df <- data.table::fread(i_file)
  
  print(paste0("nrow: ", nrow(df)))
  
  df %>% 
    dplyr::filter(generation %in% c(0, G_n - 1)) ->
    df_firstlast
  
  return(df_firstlast)
  
}

# Run below

start_time <- Sys.time()

for (i_dir in dir_pattern) {
  
  DataDir <- "./Results/"
  #DataDir <- "./Analysis_Script/Results/"
  target_dir <- dir(DataDir, recursive = TRUE)[stringr::str_detect(dir(DataDir, recursive = TRUE), pattern = i_dir)]
  All_csv <- paste0(DataDir, target_dir[stringr::str_detect(target_dir, pattern = "_nat.csv.gz")])
  #print(All_csv)
  
  binded_df <- purrr::map_dfr(All_csv, ~extract_FirstLast(.x))
  
  save_dir <- dir(DataDir)[stringr::str_detect(dir(DataDir), pattern = i_dir)]
  file_name <- "/first_lastgene.csv"
  
  readr::write_csv(x = binded_df, path = paste0(DataDir, save_dir, file_name))
  
}

end_time <- Sys.time()
print(paste0("end - start: ", difftime(end_time, start_time))) # display the elapsed time

