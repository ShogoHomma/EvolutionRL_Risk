# read simulation data and make a single data frame

MakeDf <- function(dir_pattern, sim_name, DataDir, DataDir_sub, csv_pattern) {
  
  data_path <- make_DataPath(DataDir, DataDir_sub, dir_pattern, csv_pattern)
  
  df <- purrr::map_dfr(data_path, ~data.table::fread(.x))
  
  # for RandMultiTask
  if ("sim" %in% colnames(df) ) {
    #print("found sim col !")
    df <- df %>% 
      dplyr::rename(
        sim_i = sim
      )
  }
  
  df <- df %>%
    dplyr::mutate(sim = sim_name) %>%
    dplyr::select(sim, everything()) 
  
  print(dir_pattern)
  return(df)
  
}
