# make paths for one or multiple data file (e.g., .csv)
## Dir：target directory, directly under the project directory
## Dir_sub：sub directory under target directory

#library(here)

make_DataPath <- function(Dir, Dir_sub, dir_pattern, csv_pattern) {
  
  target_path <- here(Dir, Dir_sub)
  target_dir <- dir(target_path, recursive = TRUE)[stringr::str_detect(dir(target_path, recursive = TRUE), pattern = dir_pattern)]
  target_csvs <- here(target_path, target_dir[stringr::str_detect(target_dir, pattern = csv_pattern)])
  
  return(target_csvs)
  
}
