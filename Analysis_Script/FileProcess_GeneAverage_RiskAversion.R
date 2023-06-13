# Script for making GeneAverage_*.csv
# read _nat.csv.gzip, and calculate mean value and SD of rate of choosing option 1 and probability choosing op 1 for each generation

library(tidyverse)
library(data.table)

variance <- function(x) var(x)*(length(x)-1)/length(x)

# e.g., dir_pattern <- c( "_20220322151600", "_20220322151601", "_20220322151602", "_20220323153501", "_20220323153502")
dir_pattern <- c()

for (i_dir in dir_pattern) {
  
  DataDir <- "./Results/"
  target_dir <- dir(DataDir, recursive = TRUE)[stringr::str_detect(dir(DataDir, recursive = TRUE), pattern = i_dir)]
  All_csv <- paste0(DataDir, target_dir[stringr::str_detect(target_dir, pattern = ".gz")])
  #print(All_csv)
  
  for (i_file in All_csv) {
    
    index <- which(All_csv == i_file)
    df <- data.table::fread(i_file)
    print(i_file)
    print(paste0("nrow: ", nrow(df)))
    
    df_summary <- 
      df %>% 
      dplyr::group_by(replication, generation) %>% # for each replication and generation
      dplyr::summarise(
        MEAN_payoff = mean(mean_payoff),
        SD_payoff = sqrt(variance(mean_payoff)),
        MEAN_p = mean(mean_p),
        SD_p = sqrt(variance(mean_p)),
        MEAN_rate_choice1 = mean(rate_choice1),
        SD_rate_choice1 = sqrt(variance(rate_choice1))
      )
    
    save_dir <- dir(DataDir)[stringr::str_detect(dir(DataDir), pattern = i_dir)]
    file_name <- paste0("/GeneAverage_r00", index, ".csv")
    
    readr::write_csv(x = df_summary, path = paste0(DataDir, save_dir, file_name))
    
    # Dropboxの同名フォルダにも保存する
    #DataDir2 <- "./Results/"
    #readr::write_csv(x = df_summary, path = paste0(DataDir2, save_dir, file_name))
    
  }
  
}

