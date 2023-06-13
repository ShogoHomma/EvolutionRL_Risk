# extract the first and last generation data from _nat (agent) data of random multiple-task simulation

MakeDfNat_Rand <- function(dir_pattern, sim_name, DataDir, DataDir_sub, csv_pattern, GeneN, AgentN, TaskN, Extract_SimN) {
  
  data_path <- make_DataPath(DataDir, DataDir_sub, dir_pattern, csv_pattern)
  
  df <- purrr::map_dfr(1:Extract_SimN, ~Read_FirstLast_Rand(.x, data_path, GeneN, AgentN, TaskN)) 
  
  df <- df %>% 
    dplyr::rename(task_i = task) %>%
    dplyr::mutate(sim = sim_name) %>% 
    dplyr::select(sim, everything())
  
  return(df)
  
}
