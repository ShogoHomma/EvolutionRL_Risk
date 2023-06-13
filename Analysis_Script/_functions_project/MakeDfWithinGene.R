# read data of effect of ap and an (data within a generation) and make a single data frame

MakeDfWithinGene <- function(dir_pattern, sim_name, DataDir, DataDir_sub, csv_pattern) {
  
  df <- MakeDf(dir_pattern, sim_name, DataDir, DataDir_sub, csv_pattern) 
  
  df_rev <- df %>%
    dplyr::mutate(
      replication = rep(1:10, each = 3e4)
    ) %>%
    dplyr::select(sim, replication, everything()) %>% 
    tidyr::separate(col = "sim", into = c("m2", "sd2", "m1", "sd1"), sep = "\\.", remove = FALSE) %>% 
    dplyr::mutate(
      m1 = as.integer(m1),
      sd1 = as.integer(sd1),
      m2 = as.integer(m2),
      sd2 = as.integer(sd2)
    )
  
  #print(dir_pattern)
  return(df_rev)
  
}